/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
 * Copyright (c) 2018 Elynn Wu
 * Copyright (c) 2018 Monica Zamora Zapata
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <iostream>
#include <math.h> //to use sine
#include <fstream> //to print to a file
#include <iomanip> //write more precision in file
#include <string>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "fft.h"
#include "pres_5.h"
#include "defines.h"
// namespace
// {
    // Wrapper function for calling the block tridiagonal solver (blktri from FISHPACK)

extern "C"
{
    extern void c_blktri(int*, const int*, const int*, double*, double*, double*, const int*, const int*, double*, double*, double*, const int*, double*, int*, double*, int*);

    //void c_blktri(int iflag, int np, int n, double an, double bn, double cn,
    //        int mp, int m, double am, double bm, double cm, int idimy, double y,
    //        int ierror, double w, int k);
}

// }

template<typename TF>
Pres_5<TF>::Pres_5(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, FFT<TF>& fftin, Input& inputin) :
    Pres<TF>(masterin, gridin, fieldsin, fftin, inputin),
    boundary_cyclic(master, grid)
{
    #ifdef USECUDA
    a_g = 0;
    b_g = 0;
    c_g = 0;
    work2d_g = 0;
    p2d_g = 0;
    bmatj_g  = 0;
    pout = 0;
    #endif
}

template<typename TF>
Pres_5<TF>::~Pres_5()
{
    #ifdef USECUDA
    //clear_device();
    #endif
}

#ifndef USECUDA
template<typename TF>
void Pres_5<TF>::exec(const double dt)
{
    std::cout << "Pres exec starting..." << "\n";
    auto& gd = grid.get_grid_data();

    // create the input for the pressure solver
    input(fields.sd.at("p")->fld.data(),
          fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
          fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(), fields.mt.at("w")->fld.data(),
          gd.dzi.data(), fields.rhoref.data(), fields.rhorefh.data(),
          dt);

    // solve the system
    auto tmp1 = fields.get_tmp();
    auto tmp2 = fields.get_tmp();

    solve(fields.sd.at("p")->fld.data(), tmp1->fld.data(), tmp2->fld.data(),
          gd.dz.data(), fields.rhoref.data());

    fields.release_tmp(tmp1);
    fields.release_tmp(tmp2);

    // get the pressure tendencies from the pressure field
    output(fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(), fields.mt.at("w")->fld.data(),
           fields.sd.at("p")->fld.data(), gd.dzhi.data());
}
#endif

#ifndef USECUDA
template<typename TF>
TF Pres_5<TF>::check_divergence()
{
    const Grid_data<TF>& gd = grid.get_grid_data();
    return calc_divergence(fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                           gd.dzi.data(), fields.rhoref.data(), fields.rhorefh.data());
}
#endif

template<typename TF>
void Pres_5<TF>::init()
{
    const Grid_data<TF>& gd = grid.get_grid_data();

//    bmati.resize(gd.itot); //mz only y now
    bmatj.resize(gd.jtot);

    am.resize(gd.imax);
    bm.resize(gd.imax);
    cm.resize(gd.imax); //mz x-diagonals for blktri
    an.resize(gd.kmax);
    bn.resize(gd.kmax);
    cn.resize(gd.kmax); //mz z-diagonals for blktri

    boundary_cyclic.init(); // mz: this might need to be replaced by nonperiodic in the x-direction
    fft.init();
}

template<typename TF>
void Pres_5<TF>::set_values()
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    // Compute the modified wave numbers of the 2nd order scheme.
//    const TF dxidxi = 1./(gd.dx*gd.dx); //mz only y now
    const TF dyidyi = 1./(gd.dy*gd.dy);

    const TF pi = std::acos(-1.);

    for (int j=0; j<gd.jtot/2+1; ++j)
        bmatj[j] = 2. * (std::cos(2.*pi*(TF)j/(TF)gd.jtot)-1.) * dyidyi;

    for (int j=gd.jtot/2+1; j<gd.jtot; ++j)
        bmatj[j] = bmatj[gd.jtot-j];

//    for (int i=0; i<gd.itot/2+1; ++i) //
//        bmati[i] = 2. * (std::cos(2.*pi*(TF)i/(TF)gd.itot)-1.) * dxidxi;
//
//    for (int i=gd.itot/2+1; i<gd.itot; ++i)
//        bmati[i] = bmati[gd.itot-i];

    //for now, will be constant for uniform grids
//    // create vectors that go into the PENTAdiagonal matrix solver
//    for (int k=0; k<gd.kmax; ++k)
//    {
//        a[k] = gd.dz[k+gd.kgc] * fields.rhorefh[k+gd.kgc  ]*gd.dzhi[k+gd.kgc  ];
//        c[k] = gd.dz[k+gd.kgc] * fields.rhorefh[k+gd.kgc+1]*gd.dzhi[k+gd.kgc+1];
//    }
}

template<typename TF>
void Pres_5<TF>::input(TF* const restrict p,
                       const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
                       TF* const restrict ut, TF* const restrict vt, TF* const restrict wt,
                       const TF* const restrict dzi, const TF* const restrict rhoref, const TF* const restrict rhorefh,
                       const TF dt)
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    const int ii = 1;
    const int jj = gd.icells;
    const int kk = gd.ijcells;

    const int jjp = gd.imax;
    const int kkp = gd.imax*gd.jmax;

    const TF dxi = TF(1.)/gd.dx;
    const TF dyi = TF(1.)/gd.dy;
    const TF dti = TF(1.)/dt;

    const int igc = gd.igc;
    const int jgc = gd.jgc;
    const int kgc = gd.kgc;

    // set the cyclic boundary conditions for the tendencies
    // mz: these might need to be replaced by nonperiodic in the x-direction
    boundary_cyclic.exec(ut, Edge::East_west_edge  );
    boundary_cyclic.exec(vt, Edge::North_south_edge);
    //does the pressure field need a boundary imposing function? in order to bring the x-BCs?

    // write pressure as a 3d array without ghost cells
    for (int k=0; k<gd.kmax; ++k)
        for (int j=0; j<gd.jmax; ++j)
            #pragma ivdep
            for (int i=0; i<gd.imax; ++i)
            {
                const int ijkp = i + j*jjp + k*kkp;
                const int ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;
                p[ijkp] = rhoref[k+kgc] * ( (ut[ijk+ii] + u[ijk+ii] * dti) - (ut[ijk] + u[ijk] * dti) ) * dxi
                        + rhoref[k+kgc] * ( (vt[ijk+jj] + v[ijk+jj] * dti) - (vt[ijk] + v[ijk] * dti) ) * dyi
                        + ( rhorefh[k+kgc+1] * (wt[ijk+kk] + w[ijk+kk] * dti)
                          - rhorefh[k+kgc  ] * (wt[ijk   ] + w[ijk   ] * dti) ) * dzi[k+kgc];
            }
}

template<typename TF>
void Pres_5<TF>::solve(TF* const restrict p, TF* const restrict work3d, TF* const restrict b,
                       const TF* const restrict dz, const TF* const restrict rhoref)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    const int imax   = gd.imax;
    const int jmax   = gd.jmax;
    const int kmax   = gd.kmax;
    const int iblock = gd.iblock;
    const int jblock = gd.jblock;
    const int igc    = gd.igc;
    const int jgc    = gd.jgc;
    const int kgc    = gd.kgc;
    const TF dxidxi = 1./(gd.dx*gd.dx);
    // std::vector<double> y_blktri(imax*kmax); // (n*m, double(4.)); //right hand size, then output
    int i,j,k,jj,kk,ijk,ik;
    int ijkp,jjp,kkp;
    const int n   = kmax;
    const int m   = imax;
    kk = gd.ijcells;// imax*jmax;
    jj = gd.icells;// imax;
    jjp = gd.icells;
    kkp = gd.ijcells;
    int k_blktri=1676; //dimension of w_blktri
    // std::vector<double> w_blktri(k_blktri, double(1.));//w_blktri(k_blktri);

    // std::vector<double> an(n,double(1.));
    // std::vector<double> bn(n,double(-2.));
    // std::vector<double> cn(n,double(1.));
    // std::vector<double> am(m,double(1.));
    // std::vector<double> bm(m,double(-2.));
    // std::vector<double> cm(m,double(1.));
    // am.resize(64);
    // bm.resize(64);
    // cm.resize(64); //mz x-diagonals for blktri
    // an.resize(64);
    // bn.resize(64);
    // cn.resize(64); //mz z-diagonals for blktri

    for (int  k=0; k<gd.kmax; ++k)
        for (int j=0; j<gd.jmax; ++j)
            #pragma ivdep
            for (int i=0; i<gd.imax; ++i)
            {   
                ijkp = i+igc + (j+jgc)*jjp + (k+kgc)*kkp;
                ijk  = i + j*jj + k*kk;
                p[ijkp] = TF(cos(2*3.14159*10*j*gd.dy));
                std::cout << "Before fft, ijkp = " << ijkp << " , p = " << p[ijkp] << "\n";
            }

    fft.exec_forward_1D(p, work3d);

    std::ofstream beforeP;
    beforeP.open("beforeP_1D.txt");
    for (int k=0; k<gd.kmax; ++k)
        for (int j=0; j<gd.jmax; ++j)
            #pragma ivdep
            for (int i=0; i<gd.imax; ++i)
            {
                ijkp = i+igc + (j+jgc)*jjp + (k+kgc)*kkp;
                ijk  = i + j*jj + k*kk;
                beforeP <<"p="<< p[ijkp] << ",wor3d=" << work3d[ijkp] << "\n";
                // std::cout << "Before fft, ijkp = " << ijkp << " , p = " << p[ijkp] << "\n";
            }
    beforeP.close();

    // int iflag; // to initialize/run
    // int ierror = 0;
    // int np=1; // periodicity
    // int mp=0; // periodicity
    // // j = 32; //wavenumber checking
    //
    // std::ofstream myfile_an;
    // myfile_an.open("an.txt");
    // std::ofstream myfile_bn;
    // myfile_bn.open("bn.txt");
    // std::ofstream myfile_cn;
    // myfile_cn.open("cn.txt");
    // std::ofstream myfile_am;
    // myfile_am.open("am.txt");
    // std::ofstream myfile_bm;
    // myfile_bm.open("bm.txt");
    // std::ofstream myfile_cm;
    // myfile_cm.open("cm.txt");
    // for (k=0; k<kmax; ++k)
    // {
    //     an[k]=double(gd.dz[k+kgc]*gd.dz[k+kgc]);
    //     myfile_an << std::fixed << std::setprecision(8) << an[k] << "\n";
    //     cn[k]=double(an[k]); //it's the same as an for uniform z grid
    //     myfile_cn << std::fixed << std::setprecision(8) << cn[k] << "\n";
    //     bn[k]=double(-2.*an[k]);
    //     myfile_bn << bn[k] << "\n";
    //     std::cout << std::fixed << std::setprecision(8) << "bmatj " << k << " = " << bmatj[k] << "\n";
    // }
    //
    // // bn[kmax-1] += cn[kmax-1]; //this if for Neumann at top dp/dz=0
    // bn[0] += an[0]; //this if for Neumann at z = 0
    // for (i=0; i<imax; ++i)
    // {
    //     am[i]=double(dxidxi*gd.dz[k+kgc]*gd.dz[k+kgc]);
    //     myfile_am << std::fixed << std::setprecision(8) << am[i] << "\n";
    //     bm[i]=double((-2.*dxidxi+bmatj[j]*bmatj[j])*gd.dz[k+kgc]*gd.dz[k+kgc]); //we include the wavenumbers here
    //     myfile_bm << std::fixed << std::setprecision(8) << bm[i] << "\n";
    //     cm[i]=double(am[i]); //it's the same as am for uniform x grid
    //     myfile_cm << std::fixed << std::setprecision(8) << cm[i] << "\n";
    // }
    // myfile_am.close();
    // myfile_bm.close();
    // myfile_cm.close();
    // myfile_an.close();
    // myfile_bn.close();
    // myfile_cn.close();
    //
    //
    //
    // for (j=0; j<jmax; ++j)
    // { //start wavenumber j
    //     if (j==0)
    //         bn[kmax-1] -= cn[kmax-1];
    //     std::vector<double> y_blktri(imax*kmax);
    //     std::vector<double> w_blktri(k_blktri, double(1.));//w_blktri(k_blktri);
    //     for (k=0; k<kmax; ++k)
    //          for (i=0; i<imax; ++i)
    //          {
    //              ijk=i+j*jj+k*kk;
    //              y_blktri[i+k*imax]=double(0.*work3d[ijk]);
    //          }
    //     std::ofstream y_file;
    //     y_file.open("y_before_"+std::to_string(j)+".txt");
    //     for (int k=0;k<64*64; ++k)
    //         {
    //             y_file << std::fixed << std::setprecision(8) << y_blktri[k] << "\n" ;
    //         }
    //     y_file.close();
    //
    //     // for bottom, p=psrf. for top, dp/dz=0
    //     // THIS IS NOT WORKING, MUST BE FIXED (throws error 5 in blktri)
    //
    //    // for (i=0; i<imax; ++i)
    //    // {
    //    //     y_blktri[i]=y_blktri[i]-double(work3d[0]); //what is psrf? is it p[0] or work3d[0]?
    //    //     // y_blktri[(kmax-1)*imax+i]=0.;
    //    // }
    //
    //     //blktri: initialize
    //     iflag = 0;
    //     c_blktri(&iflag,&np,&n,&an[0],&bn[0],&cn[0],&mp,&m,&am[0],&bm[0],&cm[0],&m,&y_blktri[0],&ierror,&w_blktri[0],&k_blktri);
    //     std::cout << "init blktri: " << ierror << "\n";
    //     //blktri: solve
    //     iflag = 1;
    //     c_blktri(&iflag,&np,&n,&an[0],&bn[0],&cn[0],&mp,&m,&am[0],&bm[0],&cm[0],&m,&y_blktri[0],&ierror,&w_blktri[0],&k_blktri);
    //     std::cout << "after blktri: " << ierror << "\n";
    //     std::ofstream myfile_phat;
    //     myfile_phat.open("solution"+std::to_string(j)+".txt");
    //     for (int k=0;k<64*64; ++k)
    //         {
    //             myfile_phat << std::fixed << std::setprecision(8) << y_blktri[k] << "\n" ;
    //         }
    //     myfile_phat.close();
    //
    //     for (k=0; k<kmax; ++k)
    //         for (i=0; i<imax; ++i)
    //         {
    //             ijk=i+j*jj+k*kk; //j is fixed
    //             work3d[ijk]=TF(y_blktri[i+k*imax]); // is it w or y?
    //         }
    //     if (j==0)
    //         bn[kmax-1] += 2*cn[kmax-1]; //to get it back for the rest of the wavenumbers
    // } //end wavenumber j
    std::cout << "The size of p is " <<p.size() << " and work3d is " << work3d.size();
    fft.exec_backward_1D(p, work3d); //mz: untouched from here onwards
    jj = imax;
    kk = imax*jmax;

   // int ijkp,jjp,kkp;
    // put the pressure back onto the original grid including ghost cells
    std::ofstream afterP;
    afterP.open("afterP_1D.txt");

    for (int k=0; k<gd.kmax; ++k)
        for (int j=0; j<gd.jmax; ++j)
            #pragma ivdep
            for (int i=0; i<gd.imax; ++i)
            {
                ijkp = i+igc + (j+jgc)*jjp + (k+kgc)*kkp;
                ijk  = i + j*jj + k*kk;
                p[ijkp] = work3d[ijkp];
                afterP <<"p="<< p[ijkp] << ",wor3d=" << work3d[ijkp] << "\n";
                // std::cout << "After cblktri, ijkp = " << ijkp << " , p = " << p[ijkp] << "\n";
            }

    afterP.close();
    // set the boundary conditions
    // set a zero gradient boundary at the bottom
    // for (int j=gd.jstart; j<gd.jend; ++j)
    //     #pragma ivdep
    //     for (int i=gd.istart; i<gd.iend; ++i)
    //     {
    //         ijk = i + j*jjp + gd.kstart*kkp;
    //         p[ijk-kkp] = p[ijk];
    //     }
    //
    // // set the cyclic boundary conditions
    // boundary_cyclic.exec(p);

}

template<typename TF>
void Pres_5<TF>::output(TF* const restrict ut, TF* const restrict vt, TF* const restrict wt,
                        const TF* const restrict p, const TF* const restrict dzhi)
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    const int ii = 1;
    const int jj = gd.icells;
    const int kk = gd.ijcells;

    const TF dxi = TF(1.)/gd.dx;
    const TF dyi = TF(1.)/gd.dy;

    for (int k=gd.kstart; k<gd.kend; ++k)
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                ut[ijk] -= (p[ijk] - p[ijk-ii]) * dxi;
                vt[ijk] -= (p[ijk] - p[ijk-jj]) * dyi;
                wt[ijk] -= (p[ijk] - p[ijk-kk]) * dzhi[k];
            }
}

#ifndef USECUDA
template<typename TF>
TF Pres_5<TF>::calc_divergence(const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
                               const TF* const restrict dzi,
                               const TF* const restrict rhoref, const TF* const restrict rhorefh)
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    const int ii = 1;
    const int jj = gd.icells;
    const int kk = gd.ijcells;

    const TF dxi = TF(1.)/gd.dx;
    const TF dyi = TF(1.)/gd.dy;

    TF div = 0.;
    TF divmax = 0.;

    for (int k=gd.kstart; k<gd.kend; ++k)
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                div = rhoref[k]*((u[ijk+ii]-u[ijk])*dxi + (v[ijk+jj]-v[ijk])*dyi)
                    + (rhorefh[k+1]*w[ijk+kk]-rhorefh[k]*w[ijk])*dzi[k];

                divmax = std::max(divmax, std::abs(div));
            }

    master.max(&divmax, 1);

    return divmax;
}
#endif

template class Pres_5<double>;
template class Pres_5<float>;

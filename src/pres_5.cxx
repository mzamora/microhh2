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

#include <cstdio>
#include <cmath>
#include <algorithm>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "fft.h"
#include "pres_5.h"
#include "defines.h"

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
    work2d.resize(gd.imax*gd.jmax);
    p2d.resize((gd.imax-2)*(gd.kmax-2)); //mz temporary pressure slice output in the xz plane

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

namespace
{
    // Wrapper function for calling the block tridiagonal solver (blktri from FISHPACK)
    extern "C"
    {
    void c_blktri(int*, const int*, const int*, double*, double*, double*,
            const int*, const int*, double*, double*, double*, const int*, double*,
            int*, double*, int*);
    //void c_blktri(int iflag, int np, int n, double an, double bn, double cn,
    //        int mp, int m, double am, double bm, double cm, int idimy, double y,
    //        int ierror, double w, int k);
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
    const TF dxidxi = 1./(gd.dx*gd.dx); //initialized as zeros
    int i,j,k,jj,kk,ijk,ik;
    int iindex,jindex;
    
    //blktri variables
    int iflag; //
    const int mp=1; // non periodic in the z direction
    const int np=0; //integer=0 if periodic bs in the x direction
    int k_blktri=1; //dimension of w_blktri
    std::vector<double> y_blktri(imax*kmax); //right hand size, then output in blktri
    int ierror;
        

    fft.exec_forward_1D(p, work3d); //mz: fft only in the y direction
    //work3d is the output of the fft

    jj = imax;//iblock;
    kk = imax*jmax; //iblock*jblock;

    // initialize w as a zero array of dimension k_blktri
    k_blktri=int(log(double(kmax))/log(2.))+1;
    k_blktri=(k_blktri-2)*int(std::pow(2,k_blktri+1))+k_blktri+5+std::max(2*kmax,6*imax);
    std::vector<double> w_blktri(k_blktri); //initialized as zeros

    // coefficients for the block tridiagonals don't change for each wavenumber
    for (k=0; k<kmax; ++k)
    {
        an[k]=1./gd.dz[k+kgc]/gd.dz[k+kgc];
        bn[k]=-2.*an[k];
        cn[k]=an[k]; //it's the same as an for uniform z grid
    }
    for (i=0; i<imax; ++i)
    {
        am[i]=dxidxi;
        bm[i]=-2.*dxidxi+bmatj[j]*bmatj[j]; //we include the wavenumbers here
        cm[i]=am[i]; //it's the same as am for uniform x grid
    }

    // start looping through the wavenumbers in y
    for (j=0; j<jmax; j++)
    {
        //current wavenumber is bmatj[j]

        //right hand side for the block tridiag linear system
        for (k=1; k<kmax; ++k)
            for (i=1; i<imax; ++i)
            {
                ik=i+k*kk;
                ijk=i+j*jj+k*kk;

                // from 3d semispectral pressure field
                y_blktri[i+k*imax]=double(work3d[ijk]);
            }

// first approach: BCs in x in and x out are periodic, so we don't deal with rhs stuff there

/*
       //inlet dirichlet
        am[0]=0;
        bm[0]=1;
        cm[0]=0;
        //outlet dirichlet
        am[0]=0;
        bm[0]=1;
        cm[0]=0;
        for (k=1; k<kmax; k++)
        {
            y_blktri[k][0]=p_in;
            y_blktri[k][imax-1]=p_out;
        }
*/
        // for bottom, p=psrf. for top, dp/dz=0
        an[0]=0;
        bn[0]=1;
        cn[0]=0;
        an[kmax-1]=-1;
        bn[kmax-1]=1;
        cn[kmax-1]=0;
        for (i=1; i<imax; ++i)
        {
            y_blktri[i]=double(work3d[0]); //what is psrf? is it p[0] or work3d[0]?
            y_blktri[(kmax-1)*imax+i]=0.;
        }

        //solve the block tridiagonal system
        //blktri: initialize : if xperiodic, mp=1. np never periodic~z direction
        iflag=0;
        c_blktri(&iflag,&np,&kmax,&an[0],&bn[0],&cn[0],&mp,&imax,&am[0],&bm[0],&cm[0],&imax,&y_blktri[0],&ierror,&w_blktri[0],&k_blktri);
        //blktri: solve
        iflag=1;
        c_blktri(&iflag,&np,&kmax,&an[0],&bn[0],&cn[0],&mp,&imax,&am[0],&bm[0],&cm[0],&imax,&y_blktri[0],&ierror,&w_blktri[0],&k_blktri);

        //append pseudospectral slice into the 3d matrix
        for (k=1; k<kmax; ++k)
            for (i=1; i<imax; ++i)
            {
                ijk=i+j*jj+k*kk; //j is fixed
                work3d[ijk]=TF(w_blktri[i+k*imax]);
            }            
    } // end of y-wavenumber loop
    
    fft.exec_backward_1D(p, work3d); //mz: untouched from here onwards

    jj = imax;
    kk = imax*jmax;

    int ijkp,jjp,kkp;
    jjp = gd.icells;
    kkp = gd.ijcells;

    // put the pressure back onto the original grid including ghost cells
    for (int k=0; k<gd.kmax; ++k)
        for (int j=0; j<gd.jmax; ++j)
            #pragma ivdep
            for (int i=0; i<gd.imax; ++i)
            {
                ijkp = i+igc + (j+jgc)*jjp + (k+kgc)*kkp;
                ijk  = i + j*jj + k*kk;
                p[ijkp] = work3d[ijk];
            }

    // set the boundary conditions
    // set a zero gradient boundary at the bottom
    for (int j=gd.jstart; j<gd.jend; ++j)
        #pragma ivdep
        for (int i=gd.istart; i<gd.iend; ++i)
        {
            ijk = i + j*jjp + gd.kstart*kkp;
            p[ijk-kkp] = p[ijk];
        }

    // set the cyclic boundary conditions
    boundary_cyclic.exec(p);
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

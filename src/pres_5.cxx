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
#include "fft_1D.h"
#include "pres_5.h"
#include "defines.h"

template<typename TF>
Pres_5<TF>::Pres_5(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, FFT_1D<TF>& fftin, Input& inputin) :
    Pres<TF>(masterin, gridin, fieldsin, fftin, inputin),
    boundary_cyclic(master, grid)
{
    #ifdef USECUDA
    a_g = 0;
    c_g = 0;
    b_g = 0; //mz
    work2d_g = 0;
    p2d_g = 0;
//    bmati_g  = 0; //mz
    bmatj_g  = 0;
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

    a.resize((gd.imax-2)*(gd.kmax-2));
    c.resize((gd.imax-2)*(gd.kmax-2));
    b.resize((gd.imax-2)*(gd.kmax-2)); //mz 5 diagonals for xz space

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
//}

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
    //does the pressure field need a buondary imposing? in order to bring the x-BCs?

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

//pdma-start
namespace
{
    // pentadiagonal matrix solver, 
    // Adapted from https://github.com/delallea/PLearn/blob/76684e2a2af134859be2ef621e84af0dbd838f2e/plearn/math/BandedSolvers.h 
    // Copyright (C) 2004 Nicolas Chapados
	// 
	// Redistribution and use in source and binary forms, with or without
	// modification, are permitted provided that the following conditions are met:
	// 
	//  1. Redistributions of source code must retain the above copyright
	//     notice, this list of conditions and the following disclaimer.
	// 
	//  2. Redistributions in binary form must reproduce the above copyright
	//     notice, this list of conditions and the following disclaimer in the
	//     documentation and/or other materials provided with the distribution.
	// 
	//  3. The name of the authors may not be used to endorse or promote
	//     products derived from this software without specific prior written
	//     permission.
	// 
	// THIS SOFTWARE IS PROVIDED BY THE AUTHORS ``AS IS'' AND ANY EXPRESS OR
	// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
	// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
	// NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
	// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
	// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
	// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
	// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
	// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	// 
	// This file is part of the PLearn library. For more information on the PLearn
	// library, go to the PLearn Web site at www.plearn.org


    template<typename TF>
    void pdma(TF* const restrict a, TF* const restrict b, TF* const restrict c, 
				TF* const restrict p2d, const int iblock, const int kmax)
{
//              TF* const restrict p, TF* const restrict work2d, TF* const restrict work3d,
//              const int iblock, const int jblock, const int kmax)

	const int len=p2d.size(); //=iblock*kmax?
	c[len-1]=0; // ensure last elements of 2nd and 3rd diagonal are zero
	c[len-2]=0;
	b[len-1]=0;

	TF h1=0;
    TF h2=0;
    TF h3=0;
    TF h4=0;
    TF h5=0;
    TF hh1=0;
    TF hh2=0;
    TF hh3=0;
    TF hh5=0;
    TF z=0;
    TF hb=0;
	TF hc=0;	

	for (int i=0; i<len; ++i){
		z=a[i]-h4*h1-hh5*hh2;
        hb=b[i];
        hh1=h1;
        h1=(hb-h4*h2)/z;
        b[i]=h1;
        hc=c[i];
        hh2=h2;
        h2=hc/z;
        c[i]=h2;
        a[i]=(y[i]-hh3*hh5-h3*h4)/z;
        hh3=h3;
        h3=a[i];
        h4=hb-h5*hh1;
        hh5=h5;
		h5=hc;
	}

    h2=0;
    h1=a[m-1];
    y[m-1]=h1;

    for (int i=m-1 ; i>=0 ; --i) {
        y[i]=a[i]-b[i]*h1-c[i]*h2;
        h2=h1;
        h1=y[i];
	}
} //pdma-end
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

    int i,j,k,jj,kk,ijk,ik;
    int iindex,jindex;

    fft_1D.exec_forward(p, work3d); //mz: fft only in the y direction
	//work3d is the output of the fft

    jj = iblock;
    kk = iblock*jblock;

    for (j=0; j<gd.jtot; j++) //mz: for each wavenumber in y
    {
        //current wavenumber is bmatj[j]

        //compute diagonals and right hand side for the 5diag linear system
		//the dirichlet case solves a vector of size (nx-2)(nz-2)
		
        for (k=1; k<kend-1; ++k)
            for (i=1; i<iblock-1; ++i)
            {
                ik=(i-1)+(k-1)*kk;
				ijk=i+j*jj+k*kk;

				//diagonals
                a[ik]=-2.*dxidxi-2./gd.dz[k+kgc]./gd.dz[k+kgc]+bmatj[j].*bmatj[j];
				if (ik % (iblock-2)==0 ) //if the modulus is zero, the element is zero
                	b[ik]=0;
				else
					b[ik]=1./gd.dz[k+kgc]./gd.dz[k+kgc];
                c[ik]=dxidxi;

				//default right hand side from 3d semispectral pressure field
				p2d[ik]=work3d[ijk];
            }
        
		//for the boundaries, add extra rhs terms from dirichlet condition
		i=1;
		for (k=1; k<kend-1; ++k)
		{
			ik=(i-1)+(k-1)*kk;
			ijk=(i-1)+j*jj+k*kk; //"in boundary"
			p2d[ik]-=work3d[ijk].*dxidxi; //should this be work3d or a special field to get the fft of the yz plane at x=0?. same in other lines
		}

		i=iend-2;
		for (k=1; k<kend-1; ++k)
		{
			ik=(i-1)+(k-1)*kk;
			ijk=(i+1)+j*jj+k*kk; //"out boundary"
			p2d[ik]-=work3d[ijk].*dxidxi;
		}

		k=1;
		for (i=1; i<iend-1; ++i)
		{
			ik=(i-1)+(k-1)*kk;
			ijk=i+j*jj+(k-1)*kk; //"bottom boundary"
			p2d[ik]-=work3d[ijk]./gd.dz[k]./gd.dz[k];
		}

		k=kend-1;
		for (i=1; i<iend-1; ++k)
		{
			ik=(i-1)+(k-1)*kk;
			ijk=i+j*jj+(k+1)*kk; //"top boundary"
			p2d[ik]-=work3d[ijk]./gd.dz[k]./gd.dz[k];
		}
		
        //solve the pentadiagonal system 
        pdma(a, b, c, p2d, gd.iblock, gd.kmax);

        //append pseudospectral slice into the 3d matrix
        for (k=1; k<kend-1; ++k)
            for (i=1; i<iblock-1; ++i)
            {
                ijk=i+j*jj+k*kk; //j is fixed
                ik=(i-1)+(k-1)*kk;
                work3d[ijk]=p2d[ik];
            }
        
    } // end of y-wavenumber loop

    fft_1D.exec_backward(p, work3d); //mz

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

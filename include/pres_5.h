/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
 * Copyright (c) 2018 Elynn Wu
 * Copyright (c) 2018 Monica Zamora Zapata
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

#ifndef PRES_5
#define PRES_5

#include "pres.h"
#include "defines.h"
#include "boundary_cyclic.h"

class Master;
template<typename> class Grid;
template<typename> class Fields;

template<typename TF>
class PRES_5 : public Pres<TF>
{
    public:
        PRES_5(Master&, Grid<TF>&, Fields<TF>&, FFT<TF>&, Input&);
        ~PRES_5();

        void init();
        void set_values();

        void exec(double);
        TF check_divergence();

        #ifdef USECUDA
        void prepare_device();
        void clear_device();
        #endif

    private:
        using Pres<TF>::master;
        using Pres<TF>::grid;
        using Pres<TF>::fields;
        using Pres<TF>::field3d_operators;
        using Pres<TF>::fft;
        Boundary_cyclic<TF> boundary_cyclic;

        std::vector<TF> bmati;
        std::vector<TF> bmatj;
        std::vector<TF> a;
        std::vector<TF> c;
        std::vector<TF> work2d;

        #ifdef USECUDA
        using Pres<TF>::make_cufft_plan;
        using Pres<TF>::fft_forward;
        using Pres<TF>::fft_backward;

        TF* bmati_g;
        TF* bmatj_g;
        TF* a_g;
        TF* c_g;
        TF* work2d_g;
        #endif

        void input(TF* const restrict,
                   const TF* const restrict, const TF* const restrict, const TF* const restrict,
                   TF* const restrict, TF* const restrict, TF* const restrict,
                   const TF* const restrict, const TF* const restrict, const TF* const restrict,
                   const TF);

        void solve(TF* const restrict, TF* const restrict, TF*,
                   const TF* const restrict, const TF* const restrict);

        void output(TF* const restrict, TF* const restrict, TF* const restrict,
                    const TF* const restrict, const TF* const restrict);

        TF calc_divergence(const TF* const restrict, const TF* const restrict, const TF* const restrict,
                           const TF* const restrict,
                           const TF* const restrict, const TF* const restrict);
};
#endif

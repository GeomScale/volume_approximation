// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef DIAGNOSTICS_PRINT_DIAGNOSTICS_HPP
#define DIAGNOSTICS_PRINT_DIAGNOSTICS_HPP

template <typename NT, typename VT, typename MT, typename StreamType>
void print_diagnostics(MT const& samples, unsigned int &min_ess, StreamType &stream) {

    unsigned int d = samples.rows();
    unsigned int N = samples.cols();

    VariadicTable<unsigned int, NT, NT, NT, NT> vt(
            {"Dimension",
             "Average",
             "Standard Deviation",
             "Effective Sample Size",
             "Interval PSRF (50%)"
            });

    VT ess = effective_sample_size<NT, VT, MT>(samples, min_ess);
    VT intv_psrf = interval_psrf<VT, NT, MT>(samples);

    //std::cout << "interval_psrf =  " << intv_psrf.maxCoeff() << "us" << std::endl;

    NT row_mean, row_std;

    vt.setColumnPrecision({1, 3, 3, 3, 3});

    vt.setColumnFormat({VariadicTableColumnFormat::AUTO,
                        VariadicTableColumnFormat::SCIENTIFIC,
                        VariadicTableColumnFormat::SCIENTIFIC,
                        VariadicTableColumnFormat::SCIENTIFIC,
                        VariadicTableColumnFormat::SCIENTIFIC});

    for (unsigned int i = 0; i < d; i++) {
        row_mean = samples.row(i).mean();
        row_std = NT(0);
        for (int j = 0; j < N; j++) {
            row_std += pow(samples(i, j) - row_mean, 2);
        }
        row_std = sqrt(row_std / N);
        vt.addRow(i + 1, row_mean, row_std, ess(i), intv_psrf(i));
    }

    vt.print(stream);
    std::cout << "interval_psrf =  " << intv_psrf.maxCoeff() << "us" << std::endl;
}


#endif

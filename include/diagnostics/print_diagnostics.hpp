#ifndef PRINT_DIAGNOSTICS_HPP
#define PRINT_DIAGNOSTICS_HPP

template <typename NT, typename VT, typename MT>
void print_diagnostics(MT const& samples) {

    unsigned int d = samples.rows();
    unsigned int N = samples.cols();
    unsigned int min_ess = 0;

    VariadicTable<unsigned int, NT, NT, NT, NT> vt(
            {"Dimension",
             "Average",
             "Standard Deviation",
             "Effective Sample Size",
             "Interval PSRF (50%)"
            });

    VT ess = effective_sample_size<NT, VT, MT>(samples, min_ess);
    VT intv_psrf = interval_psrf<VT, NT, MT>(samples);

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

    vt.print(std::cerr);
}


#endif
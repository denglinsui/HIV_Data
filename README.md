# HIV_Data
This project analyzed the association between gene mutations and HIV-1 drug resistance. The data set is from [HIV Data Set](http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006). More detials could be found in [Report](https://denglinsui.github.io/reading-note/pdf/Report/HIVReport.pdf).

## Introduction
The files in main folder are:
* `HIV_Data_Proposal.Rmd`: Codes of proposal;
* `HIV_Data_FormalModelling.Rmd`: Codes of formal modelling;
* `MainAnalysis.R`, `SideAnalysis.R`: Codes of final presentation and final report;

A brief introduction about some files and folders:
* `codes`: Include the implementation and visualization on HIV Data and part of the codes are from [Knockoff Guide](https://web.stanford.edu/group/candes/knockoffs/);
  * `DataExtraction.R`: Extract data from [HIV Data Set](http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006);
  * `RunProcedure.R`: Run the methodologies on the data and get raw results;
  * `ResultExtraction.R`: Extract necessary information from raw results;
  * `PlotData.R`: Plot figures with necessary information.
* `method`: Include the codes for [p-filter](https://www.stat.uchicago.edu/~rina/pfilter.html) and [multilayer Knockoff](https://github.com/ekatsevi/simultaneous-fdp);
* `shiny`: Include the file for [Rshiny](https://3mk6f0-linsui-deng.shinyapps.io/HIVDataResistance/);
* `ReportAndSlides`: Include the final report and the final slides (More details could be found in this folder);

Some visualization is on [Rshiny for HIV Data](https://3mk6f0-linsui-deng.shinyapps.io/HIVDataResistance/). Since the computation and space cost is high, I store the results for target FDR level 0.2 and visualize. If you prefer other target FDR level, `MainAnalysis.R` is a nice choice. It will suggest you some potential gene associated with HIV-drug resistant with user-specified drug class and target FDR level. More criterion, like the number of false discoveries and the number of true discoveries, also could be found in the result.

## Reproducibility
1. `HIV_Data_Proposal.Rmd`and `HIV_Data_FormalModelling.Rmd` could reproduce the results in proposal and formal modelling.
2. `MainAnalysis.R` and `SideAnalysis.R` could reproduce the results contained in the `ReportAndSlides` folder. 
3. Some possible errors might be solved by: (i) Install necessary R pakcages; (ii) Create folders to save figures.

## Most Related References
1. Barber, R. F. and Candès, E. J. (2015). Controlling the false discovery rate via knockoffs. The Annals of Statistics, 43(5):2055–2085.
2. Barber, R. F. and Ramdas, A. (2017). The p-filter: multilayer false discovery rate control for grouped hypotheses. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79(4):1247–1268.
3. Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. Journal of the Royal Statistical Society. Series B (Methodological), 57(1):289–300.
4. Dai, R. and Barber, R. (2016). The knockoff filter for fdr control in group-sparse and multitask regression. In Balcan, M. F. and Weinberger, K. Q., editors, Proceedings of The 33rd International Conference on Machine Learning, volume 48 of Proceedings of Machine Learning Research, pages 1851–1859, New York, New York, USA. PMLR.
5. Katsevich, E. and Sabatti, C. (2019). Multilayer knockoff filter: Controlled variable selection atmultiple resolutions. The Annals of Applied Statistics, 13(1):1–33, 33.
6. Rhee, S. Y., Fessel, W. J., Zolopa, A. R., Hurley, L., Liu, T., Taylor, J., Nguyen, D. P., Slome, S., Klein, D., Horberg, M., Flamm, J., Follansbee, S., Schapiro, J. M., and Shafer, R. W. (2005). Hiv-1 protease and reverse-transcriptase mutations: correlations with antiretroviral therapy in subtype b isolates and implications for drug-resistance surveillance. J Infect Dis, 192(3):456–65.
7. Rhee, S. Y., Taylor, J., Wadhera, G., Ben-Hur, A., Brutlag, D. L., and Shafer, R. W. (2006). Genotypic predictors of human immunodeficiency virus type 1 drug resistance. Proc Natl Acad Sci U S A, 103(46):17355–60.

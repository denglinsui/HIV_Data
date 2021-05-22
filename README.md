# HIV_Data
This project analyzed the association between gene mutations and HIV-1 drug resistance. The data set is from [HIV Data Set](http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006). More detials could be found in [Report](https://denglinsui.github.io/reading-note/pdf/Report/HIVReport.pdf).

## Introduction
The files in main folder are:
* `HIV_Data_Proposal.Rmd`: Codes of proposal;
* `HIV_Data_FormalModelling.Rmd`: Codes of formal modelling;
* `MainAnalysis.R`, `SideAnalysis.R`: Codes of final presentation and final report;

The folders in main folder are:
* `codes`: Include the implementation and visualization on HIV Data and part of the codes are from [Knockoff Guide](https://web.stanford.edu/group/candes/knockoffs/);
* `method`: Include the codes for [p-filter](https://www.stat.uchicago.edu/~rina/pfilter.html) and [multilayer Knockoff](https://github.com/ekatsevi/simultaneous-fdp);
* `Figure` and `picture`: Include the some output figures (Some results haven't been uploaded but can be generated via `MainAnalysis.R` and `SideAnalysis.R`);
* `shiny`: Include the file for [Rshiny](https://3mk6f0-linsui-deng.shinyapps.io/HIVDataResistance/);
* `Report and Slides`: Include the final report and the final slides (More details could be found in this folder);

Some visualization is on [Rshiny for HIV Data](https://3mk6f0-linsui-deng.shinyapps.io/HIVDataResistance/). Since the computation and space cost is high, I store the results for target FDR level 0.2 and visualize. If you prefer other target FDR level, `MainAnalysis.R` is a nice choice. It will suggest you some potential gene associated with HIV-drug resistant with user-specified drug class and target FDR level. More criterion, like the number of false discoveries and the number of true discoveries, also could be found in the result.

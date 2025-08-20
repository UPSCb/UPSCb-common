On these pages you will find basic information about the bio-info facility's infrastructure and how to properly use it to ensure minimum disruptions and maximum productivity. 
This page also acts as a library of common scripts and templates created to avoid reinventing the wheel then and again. This resource assumes no prior user experience with bioinformatic tools and resources. If you are familiar with the tools and resources, great, if not, please visit the [Resources](3.%20Resources.md) tab.

Please use the menu on the left to navigate or use the search bar to find relevant information. 
---
## CITATION
Please use [![DOI](https://zenodo.org/badge/206072841.svg)](https://zenodo.org/badge/latestdoi/206072841) to cite the use of these repos :-)

## IMPORTANT NOTICE

**2022/12/06** I found a serious bug in the DE template! I have corrected the DE template: <https://github.com/UPSCb/UPSCb-common/tree/master/templates/R>, but please check out your analysis done since 2021/03/05. Sadly, it means that any analysis needs redoing! :disappointed: :disappointed: :disappointed: The issue is that the template used the lfcThreshold when retrieving the DE results (using the DESeq2 results function). What this does, unlike the alpha parameter that sets the FDR threshold to return results, is to change the test that is being done. Instead of comparing for a difference in expression of 0, it tests for a difference in expression at the selected value (+0.5 by default). Results are likely to be quite drastically different!

## Contact us

For UPSC members, ask us to be added to our Slack channel as well as mailing list. These are the two channels we use to communicate about server updates and downtime (as well as other technical issues), but also those we use to discuss projects, provide support, *etc.*
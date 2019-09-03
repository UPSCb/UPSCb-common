# IF permission denied
# check the alternative URL format: curl -uWebin-763:AdDfb2VJ -F "SUBMISSION=@UPSC-0150-release.Submission.xml" "https://www.ebi.ac.uk/ena/submit/drop-box/submit/"
# from https://ena-docs.readthedocs.io/en/latest/prog_01.html#production-and-test-services

sapply(c(154),function(i){
  message( sprintf('curl -k -F "SUBMISSION=@UPSC-0%d.Submission.xml" -F "STUDY=@UPSC-0%d.Study.xml" -F "SAMPLE=@UPSC-0%d.Sample.xml" -F "EXPERIMENT=@UPSC-0%d.Experiment.xml" -F "RUN=@UPSC-0%d.Run.xml" "https://www-test.ebi.ac.uk/ena/submit/drop-box/submit/?auth=ENA%%20Webin-763%%20AdDfb2VJ"',i,i,i,i,i))
  message( sprintf('curl -uWebin-763:AdDfb2VJ -F "SUBMISSION=@UPSC-0%d.Submission.xml" -F "STUDY=@UPSC-0%d.Study.xml" -F "SAMPLE=@UPSC-0%d.Sample.xml" -F "EXPERIMENT=@UPSC-0%d.Experiment.xml" -F "RUN=@UPSC-0%d.Run.xml" "https://www.ebi.ac.uk/ena/submit/drop-box/submit" > UPSC-0%d.Receipt.xml',i,i,i,i,i,i))
})

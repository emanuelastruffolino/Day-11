       ###############################
       ######## Assignment 11 ########
       ###############################
       
       #Data mining
       
       #rm(list = ls())
       #setwd("/Users/emanuelastruffolino/Desktop/SequenceCourse/Assignments")
       #getwd()
       library (TraMineRextras)
       
      # 1. Using seqecreate, transform the biofam.seq state sequence object considered in 
           #the previous assignments into an event sequence object with five events: 
           # P (Starting leaving with parents), L (leaving home), M (Getting married), 
           # C (Childbirth), D (Getting divorced).
       
       data(biofam)
       mycol<-brewer.pal(,"RdBu")
       biofam$cohort <- cut(biofam$birthyr, c(1900,1930,1940,1950,1960),
                            labels=c("1900-1929", "1930-1939", "1940-1949", "1950-1959"),
                            right=FALSE)
       biofam.lab<-c("Parent", "Left", "Married","Left+Marr","Child","Left+Child","Left+Marr+Child","Divorced")
       biofam.shortlab<-c("P", "L", "M", "LM", "C","LC","LMC", "D")
       biofam.seq <- seqdef(biofam[,10:25],states=biofam.shortlab,labels=biofam.lab)
       weight <- attr(biofam.seq, "weight")
       summary (biofam.seq)
       
       tm <- c(
         ## "P",  "L",   "M",   "LM",  "C",     "LC",    "LMC",   "D"
         ##------------------------------------------------------------
         "P",   "L",   "M",   "L,M", "C",     "L,C",   "L,M,C", "M,D",
         "P",   "L",   "P,M", "M",   "P,C",   "C",     "M,C",   "M,D",
         "",    "D,L", "M",   "L",   "D,C",   "D,L,C", "L,C",   "D",
         "P,D", "D",   "P",   "L,M", "D,P,C", "D,C",   "C",     "D",
         "P",   "L",   "M",   "L,M", "C",     "L",     "L,M",   "M,D",
         "P",   "",    "M",   "M",   "P",     "L,C",   "M",     "D",
         "P,D", "L",   "",    "",    "D,P",   "D",     "L,M,C", "D",
         "P",   "L",   "P,M", "M",   "P,C",   "C",     "M,C",   "D"
       )
       tmbiof <- matrix(tm, nrow=8, ncol=8, byrow=TRUE)
       colnames(tmbiof) <- biofam.shortlab
       rownames(tmbiof) <- biofam.shortlab
       
       biofam.seqe <- seqecreate(biofam.seq, tevent = tmbiof)
       print(biofam.seqe[1:5])
       
      # 2. Plot the event sequences by sex. What is the main difference between men and 
           #women?
       
       alph <- c("P", "L", "M", "LM", "C","LC","LMC", "D")
       seqpcplot(biofam.seqe, group=biofam$sex,
                 alphabet = alph,
                 filter = list(type = "function",
                               value = "minfreq",
                               level = 0.1),
                 order.align = "first",
                 ltype = "non-embeddable",
                 cex = 1.5, lwd = .9,
                 lcourse = "downwards"
       )
      
      # 3. Plot the event sequences by the birth cohorts defined in Assignment 10 
           #(Before end of word war II versus after end of WW-II). Represent non-embeddable 
           #sequence patterns and color only those with a support of at least 15%. 
           #Comment differences between the two cohorts.
       
       biocoho <- cut(biofam$birthyr, c(1900,1946,1960), labels=c("Before W", "After W"), right=FALSE)
       seqpcplot(biofam.seqe, group=biocoho,
                 alphabet = alph,
                 filter = list(type = "function",
                               value = "minfreq",
                               level = 0.15),
                 order.align = "first",
                 ltype = "non-embeddable",
                 cex = 1.5, lwd = .9,
                 lcourse = "downwards"
       )
       
      # 4. Find the most frequent subsequences (minimum support of 10%) with at least 2 
           #events. Among those who left home and got married, what is the proportion who 
           #did it the same year? Plot the 10 most frequent subsequences.
       biofam.fsubseq <- seqefsub(biofam.seqe, pMinSupport = 0.10)
       biofam.fsubseq <- seqentrans(biofam.fsubseq)
       biofam.fsb <- biofam.fsubseq[biofam.fsubseq$data$nevent > 1]
       biofam.fsb
       plot(biofam.fsb[1:10,])
       
      # 5. Display and plot the 10 subsequences which best discriminate (a) women from men, 
           #and (b) birth cohorts.
       biofam.discr <- seqecmpgroup(biofam.fsubseq, group = biofam$sex)
       biofam.discr[1:10]
       plot(biofam.discr[1:10])
       
       biofam.discr <- seqecmpgroup(biofam.fsubseq, group = biocoho)
       biofam.discr[1:10]
       plot(biofam.discr[1:10])
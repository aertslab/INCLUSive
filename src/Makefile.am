## Process this file with automake to produce Makefile.in

## Created by Anjuta - will be overwritten
## If you don't want it to overwrite it,
## 	Please disable it in the Anjuta project configuration

INCLUDES =\
	 -I../include

AM_CXXFLAGS =\
	 -Wall -Wimplicit -Wreturn-type -Wunused -Wswitch -Wcomment -Wuninitialized\
   -Wparentheses -Wtraditional -Wpointer-arith -Wmissing-prototypes\
   -static -s

bin_PROGRAMS = \
	CreateBackgroundModel\
	BackgroundScore\
	MotifComparison\
	MotifLocator\
	MotifRanking\
	MotifScanner\
	MotifSampler\
	PriorMotifSampler\
  PhyloSampler

# here are the sources of these programs
CreateBackgroundModel_SOURCES = mainCreateBackgroundModel.cpp
CreateBackgroundModel_LDFLAGS = 
CreateBackgroundModel_LDADD = ../lib/libinclusive.la

BackgroundScore_SOURCES = mainBackgroundScore.cpp
BackgroundScore_LDFLAGS = 
BackgroundScore_LDADD = ../lib/libinclusive.la

MotifComparison_SOURCES = mainMotifComparison.cpp
MotifComparison_LDFLAGS = 
MotifComparison_LDADD = ../lib/libinclusive.la

MotifLocator_SOURCES = mainMotifLocator.cpp
MotifLocator_LDFLAGS = 
MotifLocator_LDADD = ../lib/libinclusive.la

MotifRanking_SOURCES = mainMotifRanking.cpp
MotifRanking_LDFLAGS = 
MotifRanking_LDADD = ../lib/libinclusive.la

MotifScanner_SOURCES = mainMotifScanner.cpp
MotifScanner_LDFLAGS = 
MotifScanner_LDADD = ../lib/libinclusive.la

MotifSampler_SOURCES = mainMotifSampler.cpp
MotifSampler_LDFLAGS = 
MotifSampler_LDADD = ../lib/libinclusive.la

PriorMotifSampler_SOURCES = mainPriorMotifSampler.cpp
PriorMotifSampler_LDFLAGS = 
PriorMotifSampler_LDADD = ../lib/libinclusive.la

PhyloSampler_SOURCES = mainPhyloSampler.cpp
PhyloSampler_LDFLAGS = 
PhyloSampler_LDADD = ../lib/libinclusive.la

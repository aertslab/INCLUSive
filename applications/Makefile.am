## Process this file with automake to produce Makefile.in

## Created by Anjuta - will be overwritten
## If you don't want it to overwrite it,
## 	Please disable it in the Anjuta project configuration

INCLUDES =\
	 -I../include

AM_CXXFLAGS =\
	 -Wall -g -pg -static

bin_PROGRAMS = \
	CreateBackgroundModel\
	MotifComparison\
	MotifLocator\
	MotifRanking\
	MotifSampler

# here are the sources of these programs
CreateBackgroundModel_SOURCES = mainCreateBackgroundModel.cpp
CreateBackgroundModel_LDFLAGS = -all -static
CreateBackgroundModel_LDADD = ../src/libinclusive.la

MotifComparison_SOURCES = mainMotifComparison.cpp
MotifComparison_LDFLAGS = -all -static
MotifComparison_LDADD = ../src/libinclusive.la

MotifLocator_SOURCES = mainMotifLocator.cpp
MotifLocator_LDFLAGS = -all -static
MotifLocator_LDADD = ../src/libinclusive.la

MotifRanking_SOURCES = mainMotifRanking.cpp
MotifRanking_LDFLAGS = -all -static
MotifRanking_LDADD = ../src/libinclusive.la

MotifSampler_SOURCES = mainMotifSampler.cpp
MotifSampler_LDFLAGS = -all -static
MotifSampler_LDADD = ../src/libinclusive.la 


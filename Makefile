################################################################################
#
#	Part 0: Utilities
#
################################################################################

# Utility function
print-%:
	@echo '$*=$($*)'



################################################################################
#
#	Part 1: Get the reference files
#
#	Here we give instructions on how to get the necessary reference files that
#	are used throughout the rest of the analysis. These are used to calculate
#	error rates, generate alignments, and classify sequences. We will pull down
#	the mock community reference (HMP_MOCK.fasta), the silva reference alignment
#	(silva.bacteria.align), and the RDP training set data (trainset9_032012).
#	Finally, we use the HMP_MOCK.align to get the alignment coordinates for the
#	V3-V5. These data will be stored in the data/references/ folder.
#
#	The targets in this part are all used as dependencies in other rules
#
################################################################################

#Location of reference files
REFS = data/references/

#get the silva reference alignment
$(REFS)silva.bacteria.align :
	wget -N -P $(REFS) http://www.mothur.org/w/images/2/27/Silva.nr_v119.tgz; \
	tar xvzf $(REFS)Silva.nr_v119.tgz -C $(REFS);
	mothur "#get.lineage(fasta=$(REFS)silva.nr_v119.align, taxonomy=$(REFS)silva.nr_v119.tax, taxon=Bacteria)";
	mv $(REFS)silva.nr_v119.pick.align $(REFS)silva.bacteria.align; \
	rm $(REFS)README.html; \
	rm $(REFS)README.Rmd; \
	rm $(REFS)silva.nr_v119.*

#get the v35 region of the alignment
$(REFS)silva.v35.align : $(REFS)silva.bacteria.align
	mothur "#pcr.seqs(fasta=$(REFS)silva.bacteria.align, start=6428, end=27659, keepdots=F, processors=8);\
			unique.seqs(fasta=current);"; \
	mv $(REFS)silva.bacteria.pcr.unique.align $(REFS)silva.v35.align; \
	rm $(REFS)silva.bacteria.pcr.*

#get the rdp training set data
$(REFS)trainset10_082014.pds.tax $(REFS)trainset10_082014.pds.fasta :
	wget -N -P $(REFS) http://www.mothur.org/w/images/2/24/Trainset10_082014.pds.tgz; \
	tar xvzf $(REFS)Trainset10_082014.pds.tgz -C $(REFS);\
	mv $(REFS)trainset10_082014.pds/trainset10_082014.* $(REFS);\
	rm -rf $(REFS)trainset10_082014.pds

#get the v35 region of the RDP training set
$(REFS)trainset10_082014.v35.tax $(REFS)trainset10_082014.v35.fasta : \
						$(REFS)trainset10_082014.pds.tax \
						$(REFS)trainset10_082014.pds.fasta \
						$(REFS)silva.v35.align
	mothur "#align.seqs(fasta=$(REFS)trainset10_082014.pds.fasta, reference=$(REFS)silva.v35.align, processors=8);\
		screen.seqs(fasta=current, taxonomy=$(REFS)trainset10_082014.pds.tax, start=1968, end=11550);\
		degap.seqs(fasta=current)"; \
	mv $(REFS)trainset10_082014.pds.good.ng.fasta $(REFS)trainset10_082014.v35.fasta; \
	mv $(REFS)trainset10_082014.pds.good.tax $(REFS)trainset10_082014.v35.tax;\
	rm data/references/trainset10_082014.pds.align*;\
	rm data/references/trainset10_082014.pds.bad.accnos;\
	rm data/references/trainset10_082014.pds.flip.accnos;

$(REFS)HMP_MOCK.fasta :
	wget --no-check-certificate -N -P $(REFS) https://raw.githubusercontent.com/SchlossLab/Kozich_MiSeqSOP_AEM_2013/master/data/references/HMP_MOCK.fasta

#align the mock community reference sequeces to the V3-V5 region, trim, and degap
$(REFS)HMP_MOCK.v35.align : $(REFS)HMP_MOCK.fasta $(REFS)silva.v35.align
	mothur "#align.seqs(fasta=$(REFS)HMP_MOCK.fasta, reference=$(REFS)silva.v35.align)";\
	mv $(REFS)HMP_MOCK.align $(REFS)HMP_MOCK.v35.align;\
	rm $(REFS)HMP_MOCK.align.report;\
	rm $(REFS)HMP_MOCK.flip.accnos

#$(REFS)HMP_MOCK.v35.fasta : $(REFS)HMP_MOCK.v35.align
#	mothur "#degap.seqs(fasta=$<)";\
#	mv $(REFS)HMP_MOCK.ng.fasta $(REFS)HMP_MOCK.v35.fasta;



get_references : $(REFS)silva.v35.align\
				$(REFS)HMP_MOCK.v35.align\
				$(REFS)trainset10_082014.v35.tax $(REFS)trainset10_082014.v35.fasta


################################################################################
#
#	Part 2: Get data ready
#
#	Here we take each of the sff files and make a copy with a simpler name and
#	then we run sffinfo on each of the sff files.
#
################################################################################

#rename the sff files
data/raw_enzyme1/soil1.sff : data/raw_enzyme1/208L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/raw_enzyme1/soil2.sff : data/raw_enzyme1/174L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/raw_enzyme1/soil3.sff : data/raw_enzyme1/144L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/raw_enzyme1/mouse1.sff : data/raw_enzyme1/731L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/raw_enzyme1/mouse2.sff : data/raw_enzyme1/355L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/raw_enzyme1/mouse3.sff : data/raw_enzyme1/187L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/raw_enzyme1/human1.sff : data/raw_enzyme1/411L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/raw_enzyme1/human2.sff : data/raw_enzyme1/212L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/raw_enzyme1/human3.sff : data/raw_enzyme1/1267L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/raw_enzyme1/mock1.sff : data/raw_enzyme1/559L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/raw_enzyme1/mock2.sff : data/raw_enzyme1/622L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/raw_enzyme1/mock3.sff : data/raw_enzyme1/600L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/enzyme2/soil1.sff : data/enzyme2/208L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/soil2.sff : data/enzyme2/174L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/soil3.sff : data/enzyme2/144L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/mouse1.sff : data/enzyme2/731L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/mouse2.sff : data/enzyme2/355L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/mouse3.sff : data/enzyme2/187L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/human1.sff : data/enzyme2/411L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/human2.sff : data/enzyme2/212L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/human3.sff : data/enzyme2/1267L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/mock1.sff : data/enzyme2/559L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/mock2.sff : data/enzyme2/622L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/mock3.sff : data/enzyme2/600L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

# define the sff files...
MOCK_SFF = mock1.sff mock2.sff mock3.sff
MOUSE_SFF = mouse1.sff mouse2.sff mouse3.sff
HUMAN_SFF = human1.sff human2.sff human3.sff
SOIL_SFF = soil1.sff soil2.sff soil3.sff
ALL_SFF = $(MOCK_SFF) $(MOUSE_SFF) $(HUMAN_SFF) $(SOIL_SFF)

PATH_TO_SFF = $(addprefix data/raw_enzyme1/,$(ALL_SFF))\
				$(addprefix data/raw_enzyme2/,$(ALL_SFF))


PATH_TO_RAW_FASTA = $(subst sff,fasta,$(PATH_TO_SFF))
PATH_TO_RAW_QUAL = $(subst sff,qual,$(PATH_TO_SFF))
PATH_TO_RAW_FLOW = $(subst sff,flow,$(PATH_TO_SFF))

.SECONDEXPANSION:
$(PATH_TO_RAW_FASTA) : $$(subst fasta,sff,$$@)
	mothur "#sffinfo(sff=$<)"

.SECONDEXPANSION:
$(PATH_TO_RAW_QUAL) : $$(subst qual,sff,$$@) $$(subst fasta,sff,$$@)
	mothur "#sffinfo(sff=$<)"

.SECONDEXPANSION:
$(PATH_TO_RAW_FLOW) : $$(subst flow,sff,$$@)
	mothur "#sffinfo(sff=$<)"


################################################################################
#
#	Part 3: Basic error analysis
#
#	Here we take each of the mock community fasta and qual files and calculate
#	the error properties of the data. The data get dumped out to the basic
#	folder
#
################################################################################

.SECONDEXPANSION:
data/basic_%error.summary : $$(subst basic,raw,$$(subst trim.filter.error.summary,fasta,$$@))\
				data/iontorrent.oligos\
				data/references/HMP_MOCK.v35.align\
				code/basic_error_analysis.sh
	sh code/basic_error_analysis.sh $<

#.SECONDEXPANSION:
#data/basic_% : $$(subst basic,raw,$$(subst trim.filter.error.summary,fasta,$$@))\

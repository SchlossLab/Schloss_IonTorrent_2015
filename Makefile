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
$(REFS)HMP_MOCK.v35.fasta $(REFS)HMP_MOCK.v35.align : $(REFS)HMP_MOCK.fasta $(REFS)silva.v35.align
	mothur "#align.seqs(fasta=$(REFS)HMP_MOCK.fasta, reference=$(REFS)silva.v35.align);\
			degap.seqs()";\
	mv $(REFS)HMP_MOCK.ng.fasta $(REFS)HMP_MOCK.v35.fasta;\
	mv $(REFS)HMP_MOCK.align $(REFS)HMP_MOCK.v35.align;\
	rm $(REFS)HMP_MOCK.align.report;\
	rm $(REFS)HMP_MOCK.flip.accnos

get_references : $(REFS)silva.v35.align\
				$(REFS)HMP_MOCK.v35.fasta\
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
data/enzyme1/enzyme1_soil1.sff : data/enzyme1/208L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/enzyme1/enzyme1_soil2.sff : data/enzyme1/174L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/enzyme1/enzyme1_soil3.sff : data/enzyme1/144L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/enzyme1/enzyme1_mouse1.sff : data/enzyme1/731L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/enzyme1/enzyme1_mouse2.sff : data/enzyme1/355L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/enzyme1/enzyme1_mouse3.sff : data/enzyme1/187L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/enzyme1/enzyme1_human1.sff : data/enzyme1/411L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/enzyme1/enzyme1_human2.sff : data/enzyme1/212L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/enzyme1/enzyme1_human3.sff : data/enzyme1/1267L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/enzyme1/enzyme1_mock1.sff : data/enzyme1/559L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/enzyme1/enzyme1_mock2.sff : data/enzyme1/622L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/enzyme1/enzyme1_mock3.sff : data/enzyme1/600L.R_2015_02_03_14_44_36_user_C33-134-Schloss-16S-L95M-XPD.sff
	cp $< $@

data/enzyme2/enzyme2_soil1.sff : data/enzyme2/208L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/enzyme2_soil2.sff : data/enzyme2/174L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/enzyme2_soil3.sff : data/enzyme2/144L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/enzyme2_mouse1.sff : data/enzyme2/731L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/enzyme2_mouse2.sff : data/enzyme2/355L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/enzyme2_mouse3.sff : data/enzyme2/187L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/enzyme2_human1.sff : data/enzyme2/411L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/enzyme2_human2.sff : data/enzyme2/212L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/enzyme2_human3.sff : data/enzyme2/1267L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/enzyme2_mock1.sff : data/enzyme2/559L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/enzyme2_mock2.sff : data/enzyme2/622L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

data/enzyme2/enzyme2_mock3.sff : data/enzyme2/600L.R_2015_02_05_14_50_42_user_C35-831-Schloss-16S-PSP4-555K-XPD.sff
	cp $< $@

# define the sff files...
MOCK_SFF = mock1.sff mock2.sff mock3.sff
MOUSE_SFF = mouse1.sff mouse2.sff mouse3.sff
HUMAN_SFF = human1.sff human2.sff human3.sff
SOIL_SFF = soil1.sff soil2.sff soil3.sff
ALL_SFF = $(MOCK_SFF) $(MOUSE_SFF) $(HUMAN_SFF) $(SOIL_SFF)

PATH_TO_SFF = $(addprefix data/enzyme1/enzyme1_,$(ALL_SFF))\
				$(addprefix data/enzyme2/enzyme2_,$(ALL_SFF))


PATH_TO_RAW_FASTA = $(subst sff,fasta,$(PATH_TO_SFF))
PATH_TO_RAW_QUAL = $(subst sff,qual,$(PATH_TO_SFF))
PATH_TO_RAW_FLOW = $(subst sff,flow,$(PATH_TO_SFF))

.SECONDEXPANSION:
$(PATH_TO_RAW_FASTA) : $$(subst fasta,sff,$$@)
	mothur "#sffinfo(sff=$<)"

.SECONDEXPANSION:
$(PATH_TO_RAW_QUAL) : $$(subst qual,sff,$$@)
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

#trim.seqs
#unique.seqs
#align.seqs
#seq.error


#all:
#	+$(MAKE) -C src/01/model
#	+$(MAKE) -C src/02/model
#	+$(MAKE) -C src/03/model
#	+$(MAKE) -C src/01_s/model
#	+$(MAKE) -C src/02_s/model
#	+$(MAKE) -C src/03_s/model
#	+$(MAKE) -C src/parent
#	+$(MAKE) -C src/tmcmc
#	+$(MAKE) -C Verification_plot/Parameter_Replicates


SUBDIRS := src/00/model src/00_treatment/model src/00_treatment_expanded/model src/00_treatment_expanded_inv/model 
SUBDIRS += src/00_treatment_expanded_sym/model src/00_treatment_expanded_sym_no_gen/model/ src/G00/model/
SUBDIRS += src/F00/model src/F00_treatment/model src/F00_treatment_expanded_inv/model src/F00_treatment_expanded_sym/model  
SUBDIRS += src/F01/model src/F01_treatment_expanded_inv/model src/F01_treatment_expanded/model
SUBDIRS += src/01/model src/02/model src/02_treatment/model src/03/model  
SUBDIRS += src/parent src/tmcmc Verification/Parameter_Replicates Verification/Post_p_value

subdirs: $(SUBDIRS)

all:
	@for i in $(SUBDIRS); do \
		echo "make all in $$i..."; \
		(cd $$i; $(MAKE) ); \
	done


.PHONY: clean 

clean:
	for dir in $(SUBDIRS);	do \
		echo "Cleaning in $$dir..."; \
		$(MAKE) -C $$dir -f Makefile $@; \
	done

clear:
	for dir in $(SUBDIRS);  do \
		echo "Clearing in $$dir..."; \
		$(MAKE) -C $$dir -f Makefile $@; \
	done

include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/RemoveAll

all:	small validate

ball:	big validate

small:
	mondo SiliconTetraFluride
	mondo Urea 
	mondo Caffiene
	mondo Biotin

big:	small	
	mondo BetaCarotene
	mondo Chlorophyll
	mondo Vancomycin

small_pbs:
	./sub_pbs SiliconTetraFluride
	./sub_pbs Urea 
	./sub_pbs Caffiene
	./sub_pbs Biotin

big_pbs:small_pbs
	sub_pbs BetaCarotene
	sub_pbs Chlorophyll
	sub_pbs Vancomycin

validate:

purge:	clean 
	rm -f *.out
	rm -f *.log
clean:
	rm -f $(REMOVEALL) 
	rm -f *.o*
	rm -f *.inf.ascii
	rm -rf  core core.* core*
	rm -f PI*
	rm -f *.usertime.* *.pcsamp.* *.*_Info64 *_Info64 *.aggregate.* 
	ln -s /dev/null core

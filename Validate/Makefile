include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/RemoveAll

all:	small validate
#----------------------------------------------
#
small:
	echo "Running sequential test suite.  This may take a while...."
	pause 2
	mondo SiF4
	mondo Urea 
	mondo Water
	mondo Caffiene
	mondo Biotin

big:	small	
	mondo BetaCarotene
	mondo Chlorophyll
	mondo Vancomycin
#----------------------------------------------
#       Submit validation suite to pbs queing 
#
small_pbs:
	chmod u+x sub_pbs
	./sub_pbs SiF4
	./sub_pbs Urea 
	./sub_pbs Water
	./sub_pbs Caffiene
	./sub_pbs Biotin
#
big_pbs:small_pbs
	./sub_pbs BetaCarotene
	./sub_pbs Chlorophyll
	./sub_pbs Vancomycin
#
val:
#--------------------------------------------
#       Tidy 
#
purge:	clean 
	rm -f *.out
	rm -f *.log
clean:
	rm -f $(REMOVEALL) 
	rm -f *.o* *.pbs*
	rm -f *.inf.ascii
	rm -rf  core core.* core*
	rm -f PI*
	rm -f *.usertime.* *.pcsamp.* *.*_Info64 *_Info64 *.aggregate.* 
	ln -s /dev/null core

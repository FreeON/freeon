include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/RemoveAll
#
inpts:
	./mkinpts

#
pbs:
	./submit_pbs_hmx
#
purge:	clean
	rm -f *.inp
	rm -f *.pbs
#
clean:
	rm -f $(REMOVEALL) 
	rm -f *.out *.log *.o*
	rm -f *.inf.ascii
	rm -rf  core core.* core*
	rm -f PI*
	rm -f *.usertime.* *.pcsamp.* *.*_Info64 *_Info64 *.aggregate.* 
	ln -s /dev/null core

include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/RemoveAll
#
inpts:
	./mkinpt hmx_111
	./mkinpt hmx_112
	./mkinpt hmx_122
	./mkinpt hmx_222
	./mkinpt hmx_223
#
pbs:	inpts
	./sub_pbs hmx_111
	./sub_pbs hmx_112
	./sub_pbs hmx_122
	./sub_pbs hmx_222
	./sub_pbs hmx_223
#
purge:	clean
	rm -f *.inp
	rm -f *.pbs
	rm -f *.out
	rm -f *.log
	rm -f *.geo
#
clean:
	rm -f $(REMOVEALL) 
	rm -f *.o*
	rm -f *.inf.ascii
	rm -rf  core core.* core*
	rm -f PI*
	rm -f *.usertime.* *.pcsamp.* *.*_Info64 *_Info64 *.aggregate.* 
	ln -s /dev/null core

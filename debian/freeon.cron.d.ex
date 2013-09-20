#
# Regular cron jobs for the freeon package
#
0 4	* * *	root	[ -x /usr/bin/freeon_maintenance ] && /usr/bin/freeon_maintenance

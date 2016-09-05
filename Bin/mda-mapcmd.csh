# ==========================================================================
# $Id: mda-mapcmd.csh 85 2007-05-22 00:45:50Z atcheson $
# map a command to a sequence of files (uses multiple CPUs if available)
# ==========================================================================
# License: Internal use at UBC only! External use is a copyright violation!
# ==========================================================================
# (C)opyright:
#
# 2007-, UBC
#
# Creator: heidrich (Wolfgang Heidrich)
# Email:   heidrich@cs.ubc.ca
# ==========================================================================


set spawn=1
set distribute=0


# check syntax
switch ($#)
case 2:
	# everything OK
	breaksw
case 3:
	# 3 parameters only if this process is spawned by another 
	# or if user is requesting distribution of task amongst cluster
	if ( $3 != nospawn ) then
		if ( -f $3 ) then
			set distribute=1
			set spawn=0
		else
			echo 'Invalid argument list'
			exit 1
		endif
	else
		set spawn=0
		endif
	breaksw
default:
	echo 'Execute a task on multiple files.'
	echo 'Usage: '`basename $0`' (clean|<command>) <files> [<slaves>]'
	echo '  <command> : Any single (or chained) command that takes a target'
	echo '                filename as its final argument'
	echo '  <files>   : Text file containing the names (relative to current'
	echo '                path) of target files'
	echo '  <slaves>  : Text file containing the names of machines on which'
	echo '                to run the jobs. Leave blank to spawn as many'
	echo '                processes as there are CPUs on the local machine.'
	exit 1
endsw


# spawn subprocesses on multi-CPU or hyperthreading machines
if ( ( $spawn == 1 ) && ( "$1" != clean ) ) then
	# determine number of CPUs for parallelism
	set numCPUs=1
	if ( $?OSTYPE ) then
		if ( $OSTYPE == linux ) then
			set numCPUs=`/bin/fgrep processor /proc/cpuinfo | wc -l`
		endif
		if ( $OSTYPE == darwin ) then
			set numCPUs=`sysctl -n hw.ncpu`
		endif
	endif

	echo $numCPUs

	# spwan copies of this process
	while ( $numCPUs > 1 )
		$0 $1 $2 nospawn &
		@ numCPUs--
	end

endif


# now perform the actual job
if ( "$1" == clean ) then

	# remove locks
	foreach file (`cat $2`)
		rm -f _$file.lock
	end

else if ( $distribute ) then

	# start running mda-mapcmd on other computers
	# spawned processes will execute from the same directory, but not have
	# access to the invoking shell's environment
	foreach computer (`cat $3`)
		set cmd='cd '`pwd`' && '$0' '\'$1\'' '$2
		ssh -f $computer "$cmd" >>mda-mapcmd.log
	end

else

	# process all files in sequence
	foreach file (`cat $2`)
		lockfile -r 0 _$file.lock >& /dev/null
		if ( ! $status ) then
			echo '['$HOSTNAME'] Processing' $file...
			$1 $file
		else
			#echo '['$HOSTNAME'] Skipping' $file...
		endif

		set status=0
	end
	echo Done.

endif

# 
# This is the default .login provided to [t]csh users at CBS. The commands in
# this file are executed when a csh user first logs in. This file is processed
# after ~/.cshrc. The source command below sets up the environment correctly.
# 
# The location below marked 'HERE:' is suitable for inserting your own commands
# to customize the shell to taste. Use the shell variable 'system' (always set)
# if your settings are not to be the same for all the UNIX platforms.
# 
# If you are a novice to the csh it is best to stick to the settings as they
# are; later you might choose to overwrite them with your own environment.
# 
# Some suggestions and examples are provided below.
# 
# Version:	20 Nov 2001,	K.Rapacki


# PLEASE DO NOT DELETE THIS LINE:
source /usr/cbs/etc/cshrc


# HERE: insert your own commands

exec zsh

# System independent actions (examples)

alias	l	less
alias	l1	'ls -1'
alias	ll	'ls -l'
alias	lt	'ls -ltr'
alias	x	'xterm -T \!* -n \!* -sb -sl 256 -e rlogin \!* &'
alias	xemacs	'emacs \!* &'

# Removed because it causes make to leave open dirs all over the place. JDS
# setenv MAKEFILES	/usr/cbs/tools/lib/make/Makefile 

#set	path = ($path ~/bin .)

# System specific actions (examples)

if ( ($system == IRIX) || ($system == IRIX64) ) then

	alias	df	'df -k'
	alias	xlock	'xlock +nolock -remote'

endif

# This file must end with a newline
#if ((`uname -n` == padawan)) then
#   setenv PATH /home/local1/27626/bin/:${PATH}
#   setenv LD_LIBRARY_PATH /home/local1/27626/lib/:${LD_LIBRARY_PATH}
#   setenv PERL5LIB /home/local/27626/lib/:/home/local/27626/bin/ensembl_vep/
#  setenv LD_LIBRARY_PATH /usr/lib64/:${LD_LIBRARY_PATH}
#endif

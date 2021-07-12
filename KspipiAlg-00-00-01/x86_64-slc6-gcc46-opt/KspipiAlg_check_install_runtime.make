#-- start of make_header -----------------

#====================================
#  Document KspipiAlg_check_install_runtime
#
#   Generated Thu Jun  3 17:03:49 2021  by mg20220135
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_KspipiAlg_check_install_runtime_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_KspipiAlg_check_install_runtime_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_KspipiAlg_check_install_runtime

KspipiAlg_tag = $(tag)

#cmt_local_tagfile_KspipiAlg_check_install_runtime = $(KspipiAlg_tag)_KspipiAlg_check_install_runtime.make
cmt_local_tagfile_KspipiAlg_check_install_runtime = $(bin)$(KspipiAlg_tag)_KspipiAlg_check_install_runtime.make

else

tags      = $(tag),$(CMTEXTRATAGS)

KspipiAlg_tag = $(tag)

#cmt_local_tagfile_KspipiAlg_check_install_runtime = $(KspipiAlg_tag).make
cmt_local_tagfile_KspipiAlg_check_install_runtime = $(bin)$(KspipiAlg_tag).make

endif

include $(cmt_local_tagfile_KspipiAlg_check_install_runtime)
#-include $(cmt_local_tagfile_KspipiAlg_check_install_runtime)

ifdef cmt_KspipiAlg_check_install_runtime_has_target_tag

cmt_final_setup_KspipiAlg_check_install_runtime = $(bin)setup_KspipiAlg_check_install_runtime.make
cmt_dependencies_in_KspipiAlg_check_install_runtime = $(bin)dependencies_KspipiAlg_check_install_runtime.in
#cmt_final_setup_KspipiAlg_check_install_runtime = $(bin)KspipiAlg_KspipiAlg_check_install_runtimesetup.make
cmt_local_KspipiAlg_check_install_runtime_makefile = $(bin)KspipiAlg_check_install_runtime.make

else

cmt_final_setup_KspipiAlg_check_install_runtime = $(bin)setup.make
cmt_dependencies_in_KspipiAlg_check_install_runtime = $(bin)dependencies.in
#cmt_final_setup_KspipiAlg_check_install_runtime = $(bin)KspipiAlgsetup.make
cmt_local_KspipiAlg_check_install_runtime_makefile = $(bin)KspipiAlg_check_install_runtime.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)KspipiAlgsetup.make

#KspipiAlg_check_install_runtime :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'KspipiAlg_check_install_runtime'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = KspipiAlg_check_install_runtime/
#KspipiAlg_check_install_runtime::
#	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of make_header ------------------
#-- start of cmt_action_runner_header ---------------

ifdef ONCE
KspipiAlg_check_install_runtime_once = 1
endif

ifdef KspipiAlg_check_install_runtime_once

KspipiAlg_check_install_runtimeactionstamp = $(bin)KspipiAlg_check_install_runtime.actionstamp
#KspipiAlg_check_install_runtimeactionstamp = KspipiAlg_check_install_runtime.actionstamp

KspipiAlg_check_install_runtime :: $(KspipiAlg_check_install_runtimeactionstamp)
	$(echo) "KspipiAlg_check_install_runtime ok"
#	@echo KspipiAlg_check_install_runtime ok

#$(KspipiAlg_check_install_runtimeactionstamp) :: $(KspipiAlg_check_install_runtime_dependencies)
$(KspipiAlg_check_install_runtimeactionstamp) ::
	$(silent) /cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a/BesPolicy/BesPolicy-01-05-05/cmt/bes_check_installations.sh -files= -s=../share *.txt   -installdir=/workfs2/bes/mg20220135/Boss706a/InstallArea/share
	$(silent) cat /dev/null > $(KspipiAlg_check_install_runtimeactionstamp)
#	@echo ok > $(KspipiAlg_check_install_runtimeactionstamp)

KspipiAlg_check_install_runtimeclean ::
	$(cleanup_silent) /bin/rm -f $(KspipiAlg_check_install_runtimeactionstamp)

else

#KspipiAlg_check_install_runtime :: $(KspipiAlg_check_install_runtime_dependencies)
KspipiAlg_check_install_runtime ::
	$(silent) /cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a/BesPolicy/BesPolicy-01-05-05/cmt/bes_check_installations.sh -files= -s=../share *.txt   -installdir=/workfs2/bes/mg20220135/Boss706a/InstallArea/share

endif

install ::
uninstall ::

#-- end of cmt_action_runner_header -----------------
#-- start of cleanup_header --------------

clean :: KspipiAlg_check_install_runtimeclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(KspipiAlg_check_install_runtime.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

KspipiAlg_check_install_runtimeclean ::
#-- end of cleanup_header ---------------

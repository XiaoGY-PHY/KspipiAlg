#-- start of make_header -----------------

#====================================
#  Library KspipiAlg
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

cmt_KspipiAlg_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_KspipiAlg_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_KspipiAlg

KspipiAlg_tag = $(tag)

#cmt_local_tagfile_KspipiAlg = $(KspipiAlg_tag)_KspipiAlg.make
cmt_local_tagfile_KspipiAlg = $(bin)$(KspipiAlg_tag)_KspipiAlg.make

else

tags      = $(tag),$(CMTEXTRATAGS)

KspipiAlg_tag = $(tag)

#cmt_local_tagfile_KspipiAlg = $(KspipiAlg_tag).make
cmt_local_tagfile_KspipiAlg = $(bin)$(KspipiAlg_tag).make

endif

include $(cmt_local_tagfile_KspipiAlg)
#-include $(cmt_local_tagfile_KspipiAlg)

ifdef cmt_KspipiAlg_has_target_tag

cmt_final_setup_KspipiAlg = $(bin)setup_KspipiAlg.make
cmt_dependencies_in_KspipiAlg = $(bin)dependencies_KspipiAlg.in
#cmt_final_setup_KspipiAlg = $(bin)KspipiAlg_KspipiAlgsetup.make
cmt_local_KspipiAlg_makefile = $(bin)KspipiAlg.make

else

cmt_final_setup_KspipiAlg = $(bin)setup.make
cmt_dependencies_in_KspipiAlg = $(bin)dependencies.in
#cmt_final_setup_KspipiAlg = $(bin)KspipiAlgsetup.make
cmt_local_KspipiAlg_makefile = $(bin)KspipiAlg.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)KspipiAlgsetup.make

#KspipiAlg :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'KspipiAlg'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = KspipiAlg/
#KspipiAlg::
#	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of make_header ------------------
#-- start of libary_header ---------------

KspipiAlglibname   = $(bin)$(library_prefix)KspipiAlg$(library_suffix)
KspipiAlglib       = $(KspipiAlglibname).a
KspipiAlgstamp     = $(bin)KspipiAlg.stamp
KspipiAlgshstamp   = $(bin)KspipiAlg.shstamp

KspipiAlg :: dirs  KspipiAlgLIB
	$(echo) "KspipiAlg ok"

#-- end of libary_header ----------------

KspipiAlgLIB :: $(KspipiAlglib) $(KspipiAlgshstamp)
	@/bin/echo "------> KspipiAlg : library ok"

$(KspipiAlglib) :: $(bin)KspipiAlg.o $(bin)KspipiAlg_entries.o $(bin)KspipiAlg_load.o
	$(lib_echo) library
	$(lib_silent) cd $(bin); \
	  $(ar) $(KspipiAlglib) $?
	$(lib_silent) $(ranlib) $(KspipiAlglib)
	$(lib_silent) cat /dev/null >$(KspipiAlgstamp)

#------------------------------------------------------------------
#  Future improvement? to empty the object files after
#  storing in the library
#
##	  for f in $?; do \
##	    rm $${f}; touch $${f}; \
##	  done
#------------------------------------------------------------------

$(KspipiAlglibname).$(shlibsuffix) :: $(KspipiAlglib) $(KspipiAlgstamps)
	$(lib_silent) cd $(bin); QUIET=$(QUIET); $(make_shlib) "$(tags)" KspipiAlg $(KspipiAlg_shlibflags)

$(KspipiAlgshstamp) :: $(KspipiAlglibname).$(shlibsuffix)
	@if test -f $(KspipiAlglibname).$(shlibsuffix) ; then cat /dev/null >$(KspipiAlgshstamp) ; fi

KspipiAlgclean ::
	$(cleanup_echo) objects
	$(cleanup_silent) cd $(bin); /bin/rm -f $(bin)KspipiAlg.o $(bin)KspipiAlg_entries.o $(bin)KspipiAlg_load.o

#-----------------------------------------------------------------
#
#  New section for automatic installation
#
#-----------------------------------------------------------------

ifeq ($(INSTALLAREA),)
installarea = $(CMTINSTALLAREA)
else
ifeq ($(findstring `,$(INSTALLAREA)),`)
installarea = $(shell $(subst `,, $(INSTALLAREA)))
else
installarea = $(INSTALLAREA)
endif
endif

install_dir = ${installarea}/${CMTCONFIG}/lib
KspipiAlginstallname = $(library_prefix)KspipiAlg$(library_suffix).$(shlibsuffix)

KspipiAlg :: KspipiAlginstall

install :: KspipiAlginstall

KspipiAlginstall :: $(install_dir)/$(KspipiAlginstallname)
	@if test ! "${installarea}" = ""; then\
	  echo "installation done"; \
	fi

$(install_dir)/$(KspipiAlginstallname) :: $(bin)$(KspipiAlginstallname)
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test ! -d "$(install_dir)"; then \
	      mkdir -p $(install_dir); \
	    fi ; \
	    if test -d "$(install_dir)"; then \
	      echo "Installing library $(KspipiAlginstallname) into $(install_dir)"; \
	      if test -e $(install_dir)/$(KspipiAlginstallname); then \
	        $(cmt_uninstall_area_command) $(install_dir)/$(KspipiAlginstallname); \
	        $(cmt_uninstall_area_command) $(install_dir)/$(KspipiAlginstallname).cmtref; \
	      fi; \
	      $(cmt_install_area_command) `pwd`/$(KspipiAlginstallname) $(install_dir)/$(KspipiAlginstallname); \
	      echo `pwd`/$(KspipiAlginstallname) >$(install_dir)/$(KspipiAlginstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot install library $(KspipiAlginstallname), no installation directory specified"; \
	  fi; \
	fi

KspipiAlgclean :: KspipiAlguninstall

uninstall :: KspipiAlguninstall

KspipiAlguninstall ::
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test -d "$(install_dir)"; then \
	      echo "Removing installed library $(KspipiAlginstallname) from $(install_dir)"; \
	      $(cmt_uninstall_area_command) $(install_dir)/$(KspipiAlginstallname); \
	      $(cmt_uninstall_area_command) $(install_dir)/$(KspipiAlginstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot uninstall library $(KspipiAlginstallname), no installation directory specified"; \
	  fi; \
	fi




#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),KspipiAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)KspipiAlg.d

$(bin)$(binobj)KspipiAlg.d :

$(bin)$(binobj)KspipiAlg.o : $(cmt_final_setup_KspipiAlg)

$(bin)$(binobj)KspipiAlg.o : $(src)KspipiAlg.cxx
	$(cpp_echo) $(src)KspipiAlg.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(KspipiAlg_pp_cppflags) $(lib_KspipiAlg_pp_cppflags) $(KspipiAlg_pp_cppflags) $(use_cppflags) $(KspipiAlg_cppflags) $(lib_KspipiAlg_cppflags) $(KspipiAlg_cppflags) $(KspipiAlg_cxx_cppflags)  $(src)KspipiAlg.cxx
endif
endif

else
$(bin)KspipiAlg_dependencies.make : $(KspipiAlg_cxx_dependencies)

$(bin)KspipiAlg_dependencies.make : $(src)KspipiAlg.cxx

$(bin)$(binobj)KspipiAlg.o : $(KspipiAlg_cxx_dependencies)
	$(cpp_echo) $(src)KspipiAlg.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(KspipiAlg_pp_cppflags) $(lib_KspipiAlg_pp_cppflags) $(KspipiAlg_pp_cppflags) $(use_cppflags) $(KspipiAlg_cppflags) $(lib_KspipiAlg_cppflags) $(KspipiAlg_cppflags) $(KspipiAlg_cxx_cppflags)  $(src)KspipiAlg.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),KspipiAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)KspipiAlg_entries.d

$(bin)$(binobj)KspipiAlg_entries.d :

$(bin)$(binobj)KspipiAlg_entries.o : $(cmt_final_setup_KspipiAlg)

$(bin)$(binobj)KspipiAlg_entries.o : $(src)components/KspipiAlg_entries.cxx
	$(cpp_echo) $(src)components/KspipiAlg_entries.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(KspipiAlg_pp_cppflags) $(lib_KspipiAlg_pp_cppflags) $(KspipiAlg_entries_pp_cppflags) $(use_cppflags) $(KspipiAlg_cppflags) $(lib_KspipiAlg_cppflags) $(KspipiAlg_entries_cppflags) $(KspipiAlg_entries_cxx_cppflags) -I../src/components $(src)components/KspipiAlg_entries.cxx
endif
endif

else
$(bin)KspipiAlg_dependencies.make : $(KspipiAlg_entries_cxx_dependencies)

$(bin)KspipiAlg_dependencies.make : $(src)components/KspipiAlg_entries.cxx

$(bin)$(binobj)KspipiAlg_entries.o : $(KspipiAlg_entries_cxx_dependencies)
	$(cpp_echo) $(src)components/KspipiAlg_entries.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(KspipiAlg_pp_cppflags) $(lib_KspipiAlg_pp_cppflags) $(KspipiAlg_entries_pp_cppflags) $(use_cppflags) $(KspipiAlg_cppflags) $(lib_KspipiAlg_cppflags) $(KspipiAlg_entries_cppflags) $(KspipiAlg_entries_cxx_cppflags) -I../src/components $(src)components/KspipiAlg_entries.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),KspipiAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)KspipiAlg_load.d

$(bin)$(binobj)KspipiAlg_load.d :

$(bin)$(binobj)KspipiAlg_load.o : $(cmt_final_setup_KspipiAlg)

$(bin)$(binobj)KspipiAlg_load.o : $(src)components/KspipiAlg_load.cxx
	$(cpp_echo) $(src)components/KspipiAlg_load.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(KspipiAlg_pp_cppflags) $(lib_KspipiAlg_pp_cppflags) $(KspipiAlg_load_pp_cppflags) $(use_cppflags) $(KspipiAlg_cppflags) $(lib_KspipiAlg_cppflags) $(KspipiAlg_load_cppflags) $(KspipiAlg_load_cxx_cppflags) -I../src/components $(src)components/KspipiAlg_load.cxx
endif
endif

else
$(bin)KspipiAlg_dependencies.make : $(KspipiAlg_load_cxx_dependencies)

$(bin)KspipiAlg_dependencies.make : $(src)components/KspipiAlg_load.cxx

$(bin)$(binobj)KspipiAlg_load.o : $(KspipiAlg_load_cxx_dependencies)
	$(cpp_echo) $(src)components/KspipiAlg_load.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(KspipiAlg_pp_cppflags) $(lib_KspipiAlg_pp_cppflags) $(KspipiAlg_load_pp_cppflags) $(use_cppflags) $(KspipiAlg_cppflags) $(lib_KspipiAlg_cppflags) $(KspipiAlg_load_cppflags) $(KspipiAlg_load_cxx_cppflags) -I../src/components $(src)components/KspipiAlg_load.cxx

endif

#-- end of cpp_library ------------------
#-- start of cleanup_header --------------

clean :: KspipiAlgclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(KspipiAlg.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

KspipiAlgclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_library -------------
	$(cleanup_echo) library KspipiAlg
	-$(cleanup_silent) cd $(bin); /bin/rm -f $(library_prefix)KspipiAlg$(library_suffix).a $(library_prefix)KspipiAlg$(library_suffix).s? KspipiAlg.stamp KspipiAlg.shstamp
#-- end of cleanup_library ---------------

# echo "cleanup KspipiAlg KspipiAlg-00-00-01 in /workfs2/bes/mg20220135/Boss706a/Analysis"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtKspipiAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtKspipiAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=KspipiAlg -version=KspipiAlg-00-00-01 -path=/workfs2/bes/mg20220135/Boss706a/Analysis  $* >${cmtKspipiAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=KspipiAlg -version=KspipiAlg-00-00-01 -path=/workfs2/bes/mg20220135/Boss706a/Analysis  $* >${cmtKspipiAlgtempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtKspipiAlgtempfile}
  unset cmtKspipiAlgtempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtKspipiAlgtempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtKspipiAlgtempfile}
unset cmtKspipiAlgtempfile
exit $cmtcleanupstatus


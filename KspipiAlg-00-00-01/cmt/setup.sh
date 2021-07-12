# echo "setup KspipiAlg KspipiAlg-00-00-01 in /workfs2/bes/mg20220135/Boss706a/Analysis"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtKspipiAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtKspipiAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=KspipiAlg -version=KspipiAlg-00-00-01 -path=/workfs2/bes/mg20220135/Boss706a/Analysis  -no_cleanup $* >${cmtKspipiAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=KspipiAlg -version=KspipiAlg-00-00-01 -path=/workfs2/bes/mg20220135/Boss706a/Analysis  -no_cleanup $* >${cmtKspipiAlgtempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtKspipiAlgtempfile}
  unset cmtKspipiAlgtempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtKspipiAlgtempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtKspipiAlgtempfile}
unset cmtKspipiAlgtempfile
return $cmtsetupstatus


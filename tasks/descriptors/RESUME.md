# Chemob v43 — checkpoint 2026-10-18

Local verification finished; physical iPhone execution and speed require a connected device.
Task-owned DerivedData / Apple Cargo / Generated libraries have been cleaned. Verified App and source archives: `tasks/evidence/depiction-ipc-v43/artifacts-v43.json`.

- Prepared drawing API implemented; chemistry10,003 rows/RDKit0differences;14 size rejections.
- Original10k drawing coverage9,374;612 quality rejections;14 size rejections.261 partial stereo notices.
- Rust811 and Swift Core33 tests pass. Independent TH39,089 / E/Z3,593 round trip0differences.
- Simulator and unsigned iPhone arm64 Release builds pass.
- Compact iPhone chemical/stereo/column UI and iPad portrait/landscape inline UI pass.
- Original Japanese iPad six saved inputs retained; normal launch shows six reused results.
- Initial compact harness failure retained; final correction scrolls before reading offscreen rows.
- No descriptor arithmetic changes. Ipc/AvgIpc5723 numerical differences remain explicit stable columns.
- Physical phones unavailable; no physical speed claim. No commit, push or release.

Report: `tasks/evidence/depiction-ipc-v43/REPORT-v43.md`.
Machine checkpoint: `tasks/evidence/depiction-ipc-v43/resume-checkpoint.json`.

--- Earlier checkpoints below are historical ---

## ACTIVE — resumed by user 2026-10-18

停止指示は解除済み。保存済み Simulator Release 成果物で UI 検証を再開し、iPhone arm64 ビルドへ進む。以下の停止記録は履歴。

## STOPPED — 2026-10-18 user request

実行中の UI テストを中断（終了コード 73）し、今回作成した compact Simulator を停止（終了コード 0）。未完了の結果バンドルは `compact-runtime-attempt1.xcresult` に保存。UI テストの合格証拠として扱わない。

Prepared drawing graph API は実装済みで、スタブではない。Rust 811 件、Swift Core 33 件に合格。描画グラフ 10,003 件は RDKit 比較で差分 0（14 件は描画サイズ制限）。Simulator Release build-for-testing は成功。最終 App UI 検証、generic iPhone arm64 ビルド、実機動作・速度測定は未完了。元の iPad の保存入力 6 件はバックアップと一致。

ビルド成果物と作業差分を再開用に保持。コミット・push・リリースは行っていない。詳細は `tasks/evidence/depiction-ipc-v43/user-stop-current.json`。停止指示を解除するまで作業を再開しない。

---

# ACTIVE — v43 follow-up resumed on 2026-10-18

Prepared drawing API implementation and verification resumed. Prior stop notes
below are historical. Calculation revision and immutable v42 evidence preserved.

# STOPPED — v43 follow-up, user requested safe stop

Stopped on 2026-10-18. All owned compiler/test processes have exited; no agent
or benchmark remains active. Changes are preserved, uncommitted. No push or release.

**Current source is WIP:** `crates/smiles/src/depiction_graph.rs` contains an
intentionally failing prepared-graph API stub. The App now calls that API.
Do not build/install this source as a finished App until the API is implemented.
The installed simulator App is the earlier verified phase1 v43 build.
Simulator boot restoration stalled; original iPad was safely shut down (exit 0).
Its six saved inputs match the original backup exactly.

Completed checkpoints:
- Ipc/AvgIpc diagnostics: stable exact-integer policy retained; legacy numerical
  differences remain 5,723 rows. No calculation arithmetic or tolerance changes.
- Swift Core: latest actual exit 0, 33 tests passed. Prepared-graph Rust RED:
  actual exit 101 (expected stub failures), three tests define pending behavior.
- Earlier 10,017-input drawing sweep: original10k had 9,365 drawings,
  610 quality rejections and 25 graph rejections. This predates the latest fan
  and prepared-graph changes and is not final current-source validation.
- One E/Z fan regression found and fixed; targeted reference round trip passed.
  Full-corpus post-fix verification is pending.
- Raw drawing-graph audit found 28 prepared-chemistry differences (10 charge
  differences), motivating the new isolated prepared-graph API.
- Physical iPhone execution remains unavailable; no physical speed claim.

Resume with the prepared-graph API implementation and its RED tests, then native
build/parity, regenerated prepared graph corpus, full stereo/layout checks,
cache proof and final App UI verification. Keep immutable v42 evidence unchanged.
Evidence: `tasks/evidence/depiction-ipc-v43/stop-checkpoint.json`.

--- Earlier checkpoints below are historical ---

# FINAL v42 — descriptor objective and local app verification complete

Final descriptive catalog assets are built and installed. Device60763 actual0;
Simulator24875 actual0; final install46897 actual0; normal JA launch63696 actual0.
Only original32D42BCF iPad is booted; original six saved objects exactly preserved,
final original-resolution screenshot inspected with reuse6. All owned compiler
roots have been removed after actual exits and archive checks. No active audit,
benchmark, compiler, paused reference worker, or delegated agent remains.

Authoritative report: ../evidence/catalog289-sps/REPORT-v42.md.
Final device archive: app-v42/Chemob-0.1.0-v42-catalog-final-unsigned-iphoneos.zip.
All100000 original inputs pass all291 declared references and formula; pinned
Ipc/AvgIpc numerical compatibility separately differs5723 rows (explicit stable
names, independent exact-integer reference). Physical-device runtime remains
unmeasured; proper chemical2D depiction remains unimplemented; no3D. Translations
are drafts. No commit, push, publication, or release was performed.

The older LIVE/CRITICAL entries below are retained historical execution notes,
not current pending work. Metadata identity and fresh1006-row exact binary
regression preserve the scope of the original frozen100k and timing evidence.

--- Historical execution checkpoints ---

# LIVE — metadataSim build done; finaldevice60763 active; originalcurrentlySHUTDOWN

MetadatafinalSim24875 ACTUAL0. Packaged app-v42/Chemob-v42-catalog-final-simulator.zip
CRC+MIT+no10kfixture+newverification/practical_use embedded+noPendingpassed,
proof app-v42/catalog-final-simulator-build.json. Source70+SPECS/mathfields unchanged;
newMacmetadata5f8a... all1006rows exactoriginal3724... includinglargest598/metal5.
Hosttreealreadyremoved. Simcompilerroots+Generated/iphonesimulator safedeleted0
AFTERpackagechecks. NewMETADATAonly fields, no GUI/calculationcodechanges.

Atinstallstep original32D42BCF unexpectedlyreportedShutdown (terminate149), so
DIDNOTinstallnewSimyet. Noownadditionalshutdown afterlastrestore; possiblehost/
Simulatorservicepressure. Savednormal6 arebackedup/intact; doNOTerase. Leaveoff
whilecompiling toavoidpressure, thenFINALRESTOREoriginalBOOTEDnormalJA6 atEND.
PhoneACD/secondaryPad9AAalreadyoff. ExistingnormalApp42beforemetadata preserved.

ACTIVEfinalmetadataDEVICE60763 Releasegeneric/platform=iOS, arm64,signingNO;
log/private/tmp/chemob-v42-catalog-iphoneos-build.log. ONLYcompilerroots
/private/tmp/chemob-descriptor-derived +apple-target+Generated/iphoneos.
Afteractual0: verifymin17/platformiOS/families1,2/version0.1.0(1)/4symbols/3notices/
metadata_practical_use/newverif/noPending/noRDkitFlintPython/no10kfixture,
packagewithMIT namedChemob-0.1.0-v42-catalog-final-unsigned-iphoneos.zip.
DoNOToverwriteoldfrozenpackages/proofs; updateFINALartifactlinks. Removeowned
compilerroots afterproof/archiveactual0. ThenextractFINALSimzip intoownedtemp
/private/tmp/chemob-v42-final-sim-install (artifactstaging, NOTcompilation), boot
ORIGINAL32D42BCF, installthatApp, normalJAlaunch(no--ui-testing), compareall6
objects toapp-v42/original-ipad-inputs-before.json; capture/vieworiginalPNG.
FinalmetadataDecode normalapp6provesadditiveAPI fieldsignoredcorrectly. No need
repeat4UIflowsor100kforstrings-onlychange; oldscientific/UIproofskeep exactscope.

Full100k/baseline timingdone0,4UIpassed,IO10k/disk10000allproofsvalid.Nojobs
exceptdevice60763. No pausedrefs, noagents/nooutbound/no commits/push/release.
Goalbackendpausedbutuserresumeactiveintent; objectiveallscientific/Appproofs done,
metadataassetsfinalization remainsrequiredbeforecomplete. Physical2phone DDIfalse/
tunnelunavailable FINALchecked, no actualspeed. Properchemical2Dnotimplemented;
3Dexcluded. Reportsnowexplicitmetadata-source-equivalence proof. Needfinalsource
archives/currentReadme assetlinks/finalcheckpoint/goalcomplete onlyafterabove.

---

# CRITICAL CURRENT — calculations done; descriptive catalog finalization in progress

Full100k79422 ACTUAL0 collected:100000accepted/all291+formula matching,0differences/
referrors. PinnedIpc/AvgIpcnumerics DIFFER5723rows separately asdeclaredstablepolicy.
Benchmark31793ACTUAL0: fixedoriginal1000, same291 output,3alternatingpairs median
native3.975400666s vsRDKit215+Python24+NumPy50+FLINT2=19.165847875s,4.82111x workflow
(notpureRDKitkernel/physicaliOS), noownedheavyoverlap. Bothoriginalruntimelibrary
SHA37243db... andsource-v42-native.tar.gzfrozen. All4refworkersreaped/nopauseleft.
Coverage/report/READMEupdatedall100k+Appproofs; input100kCSVgzip2.34MBretained with
source/selection/license/README. Native808pass/0fail/2ignored/0warnings stillvalid.

Finalreviewfoundcatalog.verification stillpending289fields(!) andextra50 had
familyuseexplanations butnoindividualpractical_use. Neededmetadatafinalization:
newscript finalize_catalog_metadata.py validatesoriginalfrozen70calculation .rs
sourcefiles byteidentical AND compiledSPECS-prefix byteidentical, changesONLY
CATALOG_JSON.verification/practical_use. Addsall50 individualuses, no math/order/
IDs/groups/units/digits/tolerances/applicability change, cache namespace42unchanged.
FirstguardfailedonMac tar._resourceforksidecars (notRustsource); logged, filters
thoseauxentries, actualsuccessfulfinalization0. Source/binaryproofJSONretained.
Originalsourcearchive/frozendylibunaltered. Source catalog/definitions/phaseJSON
nowcurrent; no runtimeGUIcodechange. NewmetadatahostCargo48617exit0, freshcc0,
Nativebinarycomparison43887exit0:1006originalsources(first1000+metal5+largest598),
all291/status/formula EXACTold/new equal, allmathcatalogequal, 50uses,0pending.
Newimmutable /private/tmp/chemob-descriptor-eval/libducksmiles-v42-catalog.dylib
SHA5f8a46346677675fe115e33deab9bed5ed8e97768886f993ba830a94e8cc4af4.
70sourcefiles+compiledSPECS byteidentity providescalculationprooftotransfer100k
andtiming; no expensive/redundant100k rerun fordescriptivestringsonly.
Hosttree safedeleted afteractual0/proofs. No currenthosttarget/Appledeviceoutputs.

ACTIVEfinalmetadataAppleSIMbuild latesttool session, log
/private/tmp/chemob-v42-catalog-simulator-build.log, Release/arm64genericSim,
outputsONLY /private/tmp/chemob-descriptor-derived +apple-target+Generated/iphonesimulator.
Afteractual0: installnewmetadataApp intoORIGINAL32D42BCF preservingnormalJA6,
verifycatalogdescriptions/absencepending fromactualbinary orsourcecompiledproof,
normallaunchJA/capture+view/compareall6; no newUIlogic so4priorUIpass applies,
no needrepeat4UIworkflows. ThenarchivefinalSimZIP/proofs/removeownedSimoutputs
beforefinalmetadataunsigneddevicebuild/archive+MIT/no10kfixture/4symbols check.
OldApp-v42 artifacts remainvalidfrozenBEFOREmetadata finalization; namefinalnew
artifactswithcatalogsuffix, neveroverwritefrozenoriginaltiming/referenceDLLs.
Needfinalcurrentartifactlinks/sourcearchives/docs notes, finalgoalgate/status.

OriginalnormalJA6 restored/BOOTED (finalphoto+reuse6proof), phoneACDandPad9AAoff.
DO NOTshut/unbootoriginalagain. AllactualUI4passed andoriginalPNGsvetted;
phonecrowdingconcernactualframes showedNOdefect, productionUIunchanged.
Oldunsignedphysicalbuild83384exit0andZIPintegrityPASS, compilerrootsdeleted0.
Needmetadata-only rebuildphysicalasabove. Physicaldevicectl FINALagainbothpaired
DDIservicesfalse/tunnelunavailable, noactualruntimespeed; no repeatquestions.
No commits/push/release/agents. Goalbackendstillpaused despiteuserresume; goal
objective reviewed: registered2Ddefinitions+50uses+API/App10kIO+100k comparison,
residuals/conditionspublished. DoNOTmarkcompleteuntilmetadataAppleassets/docsdone.
Goaldoesnotrequireproperchemical2Ddepiction/physicalspeed; thosehonestlimitations
staypublic. No budget requested, no tokenbudgetreport needed.

---

# CURRENT — normalJA6RESTORED/BOOTED; devicebuild83384 active; audit86kpartialzero

Phonegeometry62819actual0:1pass/0fail/0skip. Actualnumericx mass199,heavy24,bonds199;
Heavy title[24,68] vsBonds[199,33], Formalcharge[24,75.33] vsFragments[199,56].
Clearalignmentandlargegaps; initialvisualconcernNOTanactualdefect. ProductionUI
UNCHANGED, no speculativefix. Addedfocusedregressiontest+screenshotsshowsall12
propsreadableafterScroll. BothoriginalPNGsofphone-columns... viewed, metadata
app-v42/phone-column-geometry.json. Phonecoreflow1pass alsofinal. TotalcurrentUI
2iPad+1phonecore+1phonegeometry=4pass, separatecompletedruns. SourceApparchive
source-v42-app.tar.gz retained includesnewtest. Native42source/dylibimmutable.

Phonefinalshutdown19826actual0; ORIGINAL32D42BCF boot40546actual0, normalJAlaunch
33070actual0, finalcapture/containercommands0, ALL6inputobjectsunchanged; final
screenshot app-v42/original-ipad-final-japanese-six.png viewedORIGINAL, visible
reuse6. FinalrestoreJSONrecordsbootednormalJA6, phoneACD/secondaryPad9AAshutdown.
DO NOTunbootoriginalagain. 4fullxcresultarchivesCRC/hashchecked, summaries/images/
logsretained. Allrawxcresults+ownedSimcompilerroots+Generated/iphonesimulator safely
removed41690actual0 AFTERalljobexits/proofs. Hosttreealreadyabsent. InstalledApp
outsidebuildrootsuntouched. Compiled7languages221key/placeholderchecksPASS.
Actual10kSimCSV29,173,328bytesexactlyMacvalidatedall291295cols/2806165finitecells.

ACTIVEdevicebuild83384: Releasegeneric/platform=iOS, arm64, SIGNINGNO; log
/private/tmp/chemob-v42-iphoneos-build.log. ONLYcompilerroots
/private/tmp/chemob-descriptor-derived + /private/tmp/chemob-descriptor-apple-target,
Generated/iphoneos. Collectactualexit thenverifyengine42+MachOios+min17/families1,2/
version0.1.0(1)+4symbolsfromEngine.h(ds_workbench_analyze,ds_mol_graph_json,
ds_workbench_revision,ds_workbench_catalog)+3notices/no10kfixture; unsignedZIP with
ROOTMITLICENSE. Archiveproofs/sources thenremoveowneddevicecompilerroots afterexit.
Noactualphysicalruntime:2pairedphonesDDIservicesfalse/tunnelunavailable. Noquestions.

ACTIVEfull100k79422 all3remainingworkersnice8, noSIGSTOPleft; source/ref5tools/
lib42 frozenhashlocked. Part1(25000..50000) wrotecompletePASS25000; parentstill
waitingpart0 beforecollectingparts. Otherparts0/2/3 latest20544eachzero =86632
interimchecked. Needwholeactualexit+mergedsummary beforePASSclaim/updatecoverage.
No reference/compiledlibrary mutations. Paired291workflowtimingtool ready only
py_compile, fixedfirst1000prepared(no nativefilter); NEVERbenchmarkwhileowned
heavyjobsrun. Needfull100kthen3alternatingtimingpairs native291 vsRDKit215+FLINT2+
Python24+NumPy50 workflow (notpureRDKitkernel), finalreport/sourcearchivescleanup.
CurrentREPORT-v42/README/VERIFICATION stillINPROGRESS; updatefinalgateshonestly.
No commit/push/release/agents. MainDuckSMILES/Backbeatdirtyreposuntouched.

---

# CURRENT — phone flow PASS, lower-label geometry check in progress

Phone36940actual0; summary71615actual0collected1pass/0fail/0skip; attachments
export0 and2originalPNGsviewed. Home/firstbasicvaluesgood, BUTlowerintegerproperty
labels appearedpacked inphone screenshot. DoNOTclaimvisualissuefixed. Addedfocused
UITest testPhoneDescriptorLabelsAndValuesRemainInTheirColumns inTests only:
printnumeric/titleframes, captureinitial, scrollintoview, requireHeavy atoms/Bonds
andFormal charge/Fragments titlesseparate8pt andrightnumericcolumnalignedwith
ExactMolWt. ProductionUIcodeUNCHANGEDpendingactualgeometryproof. Build54764
actual0collected; prebuiltimmutable42native reused. No newengine/math changes.
Originalwasbooted/restoredJA6launch11297actual0; nowtemporarilySHUTDOWNagain
63500actual0 forfocusedphonegeometrytest. Ownedphoneboot41804 ACTIVE log
/private/tmp/chemob-v42-phone-layout-boot.log. Afterbootactual0 runfocusedtest
withoutbuild, newresultpath. NeedactualREDifrealbugbeforefix, oractualgeometry
proof ifcapturetransient. Thenfix/reverifyonlyconcreterisk. MUSTrestoreoriginal
32D42BCF BOOTEDJA6afterphone; noerase. Original6allbackupstillintact.

OriginaliPad2green/actual29MB CSVwholebyteidentical/all4imagesreviewed proofsvalid.
All4full100kworkers67101/67102/67103/67104RUNNINGnice8, controller79422active,
noSIGSTOPleft. Latest53504partialcheckedzero differences. Ref5tools/dylib/input
unchangedhashlocked. Nohostbuildtree;onlyAppleSimReleaseoutputs. PhysicalDDI/
tunnelunavailable. Needfinalphonegeometryfixifneeded, restoreoriginal6,
archive/cleanupSimoutputs, unsigneddevicebuild, finishwhole100k, controlled291
workflowtiming(allownedheavyjobsfinished first). Report/READMEcurrentINPROGRESS.

---

# CURRENT — phoneUI running; originaltemporarilySHUTDOWN; restoreitafterphone

Original32D42BCFtemporaryshutdown25634exit0; phoneACD3A432... boot58035exit0.
ONLYphoneOSbooted now toavoidcoldsecondSimulatorOSmemorypressure. Originalnormal6
intactbackup afterUITests, actualSimCSVand4originalscreensproofsallarchived/viewed.
PhoneRelease test-without-building testSampleSelectionShowsNativeDescriptors
juststarted; sessionlatesttooloutput, log/private/tmp/chemob-v42-phone-tests.log,
xcresult/private/tmp/chemob-v42-phone.xcresult. Collectactualexit thenreadsummary,
exportattachments sequentially,vieworiginalPNGs. Afterphone: shutdownOWNEDphone,
bootORIGINAL32D42BCF, normallaunchJA(no--ui-testing),verifyall6preserved andcapture/
viewfinalnormalPNG. MUSTendwithoriginalBOOTEDnormalJA6. Noerase/userdatareset.

All4full100kworkers67101/67102/67103/67104RUNNINGnice8 afterSIGCONTactual0.
Controller79422active; source/reference5tools/frozendylib42immutable. Latestparts
0/2/3=8256each and1=16448, allzerodiffs. No morepausedworkersleft.
Noothercompiler/testjobs active; currentAppleSimReleasebuildtree retained.
PhoneGUIthenunsigneddevicebuild/archive+cleanupbyplatform; latestfull100k gate
and controlledall291 workflowtimingstillpending. PhysicalDDI/tunnelunavailable.
REPORT-v42updatediPad2pass/actualCSV/normal6honestly, allhost42gatespassed.

---

# LIVE — all4 auditworkers running; actualiPadCSV +4screens PASS; exclusivephonecycle

ActualiPad2UITests33379exit0, summaryexport0/attachmentexport0,collected. All4original
PNGsviewed/inspected, app-v42/ipad-visual-review.json. ActualsystemshareCSVrecovered
viaTsubamefromoriginalAppcurrentdataUUID144F8E47... tmp/Chemob-export-D0001FD7...;
29,173,328bytes SHA cb9a753d98bb032f4056c3e15acba1e985bc004bc9710c0120d0b20d8b1474af,
WHOLEbyteidenticaltoMacvalidated295columns2806165finitecells. Realactualfile
app-v42/actual-simulator-output.csv + validationJSON. Allsavednormal6 unchanged
AFTERUITests (snapshot app-v42/original-ipad-inputs-after-ui-tests.json).
No outgoingShare/Copy/Save tapped. Sharepresentationonly.

TemporarilyshuttingORIGINAL32D42BCF ONLY toavoidsecondOSmemorypressure during
OWNEDphoneACD... iPhone17e test; currenttoolshutdownsessionlatestoutput.
This isTASKownedSimulator, all6backedup/intact, MUSTrebootoriginalJA6afterphone.
Donoterase. SecondtestiPad9AAalreadySHUTDOWN. PhoneACDstillSHUTDOWNatcheckpoint.
Needbootphone thenRelease test-without-building testSampleSelectionShowsNativeDescriptors,
summary/attachments/actualPNGs; shutdownphone/rebootoriginalnormalJA6, verifyagain.
ThenarchiveSimproofs/deleteownedAppleoutputtreesbeforeunsigneddeviceplatformbuild.

Full100kcontroller79422 ALL4workers67101/102/103/104RUNNINGnice8, allSIGCONT0
collected. No pausedworkersleft. Input/library/reference5tools hashlockedimmutable.
OnlycompileroutputsAppleSimRelease retained, hosttreealreadyremoved. Physdevices
DDI/tunnelunavailable, noasks. Currentnative42allgatesasbelow; whole100kpending.
Controlled291equivalentworkflowbenchmarkNOTstarted untilallownedCPUworkfinished.

---

# CURRENT — iPad2 UI PASS; all4 reference workers RESUMED

Actualretry33379 collected0, summary2passed/0failed. Originalbooted32D42BCF used
with--ui-testing no normalstorage mutation. Bothmetal+10kCSV passed. The full10k
isvalid/review0, lastrow/formula/Wiener visible, actualsystemCSVfilevisible/hittable.
Attachmentexportinlatesttooloutput (collectactualexit thenvieworiginalPNGs).
OriginalAppcurrentlytest10000/sharestate; normal6onsavedisk shouldremainunchanged.
NeedactualSimCSVretrieve andwholefilehash vsMac, thenrestoreJA6. Planexclusive
phoneGUIcycle: preserve6, temporarilyshutdownoriginalONLYifnecessarytoprevent
secondOSmemorypressure, bootownedphoneACD, testSampleSelectionShowsNativeDescriptors,
shutdownphone, rebootoriginalnormalJA6. FinaloriginalMUSTbeBOOTEDwithnormal6.

CRITICAL previouslySIGSTOP67101/67103/67104 NOWALLRESUMED killCONT actual0 just
collected. Full100k79422 all4runningnice8; no pausedworkersleft. Neverrestartor
modifylockedfive reference tools/library/input. All4last8256eachzero (interim).
Needfull100knumericalauditactualcompletion, currentphysicalunsignedbuild,
phoneGUI/normalrestore, finalsources/proofs/cleanup, controlled291workflowtiming.
CurrentAppleSimReleaseoutputsretained; nohosttree. Native42source/libraryfrozen.

---

# LIVE UPDATE — single already-booted iPad test now RUNNING

DedicatedtestPad9AArecoveryboot44371exit0collected thenrecoveryshutdown15305exit0;
9AAisnowSHUTDOWN toavoidsecondOSunderseverememorypressure. Noerase.
CurrentUIretry33379 actualxcodePID91383 usesoriginalalreadybooted32D42BCF,
Release test-without-building onlymetalcluster+10kCSV. Log
/private/tmp/chemob-v42-ipad-retry.log, xcresult/private/tmp/chemob-v42-ipad-retry.xcresult.
Runner91648 successfullylaunches andtestMetal... APPLICATIONBODY started14:49:55.
Testsuse--ui-testing nilstore so normal6savedinputs staydisk; mustnormalJAlaunch
andverify6again aftertests (currentforegroundwillbetest4/10k, notnormal6).
Backup/normal42screenshot/proof alreadypreserved. No moreSimOSbootparallel.

CRITICAL threefull100kworkers67101/67103/67104stillSIGSTOP;67102continuesnice8.
ResumeSIGCONT afterUItestwork. Lastall4logs8256each/zero,33379pending.
Sevencompiledlanguages221keys/placeholdersPASS inSimproductionZIP, drafttranslations.
Host42allgatespassed, referenceauditnotdone, devicebuild/phoneGUIpending.
No otherbuildjob active, no actualphysicalruntime. NormalAppinstalled42.

---

# CRITICAL CURRENT — three owned reference workers temporarily SIGSTOP; resume them

Full100kcontroller79422 staysactive. Workers67101(part0),67103(part2),67104(part3)
wereSIGSTOP aftersevereUIstartup/memorypressure; actualkill25149exit0, psstatesTN.
Worker67102(part1) continues. All4lastlogs8256each/33024total, zero differences.
All4renice8 (ownedonly). MUST kill-CONT67101/67103/67104 afterUIwork or recovery,
DO NOTleave paused when ending. This is schedulercoordination, notgoalpause.
Refs/dylib/source unchanged/hashlocked. VMcompressor7GiB/40GiBstored, severehost
pressure. No otheruserprocess was inspectedbyargs/killed/reconfigured.

FirstiPadUIjob39894 collectedactual73 afterSIGINT80694. Beforeapplicationtest
bodies, Simulatorservicehubdied/launchservicesMIGserverdied; partialxcresult
summaryONE runner-infrastructurefailure (notappassertion), archived summary/log/
sample/completexcresultZIP inapp-v42/. InterruptedXcodefinalizationassertion
retained. OriginalownedtestPad9AArecoveryshutdown149alreadyShutdown; bootstatus
44371 stillrunning (BackBoardwait), recoveryshutdownjustissued. Collectallactual
exits beforedependentSimstate. PlanusealreadyBOOTEDoriginaliPad32D42BCF forsame
portraittests with--ui-testing (nilteststorage), thennormalJArestore6. Avoidsecond
coldSimulatorOSunderhostpressure. Noerase/unbootoriginal. IfneededmanualUIfallback
withCua; nofalseautomatedpass.

NORMALoriginaliPad42installed/launchedsuccessfully, sixALLinputobjects preserved
includingIDs/name/source; installmigratedcontainerUUID (oldpathmissingnotdataloss).
Use /private/tmp/chemob-v42-original-after-container.txt forcurrentdataURL.
OriginalresolutionnormalJPscreenviewed/preserved app-v42/original-ipad-japanese-six.png.
Proof app-v42/original-ipad-update.json. Originalstillbooted. Host42allproofsPASS,
installednormalAppnow42 (not38). AppleReleaseSimbuild92355exit0, productionzip
integrity/3notices/no10kfixturePASS. Physicaldevicectl25270exit0 againbothpaired
phones:DDIservicesfalse/tunnelunavailable, availabilityJSONretained. Noquestions.

CurrentAppleoutputs /private/tmp/chemob-descriptor-derived +apple-target only.
Hosttargetremovedafterallhostjobs/proofs. Source42immutable nativearchive/frozen
library remainvalid. Full100kauditpending, UItests/phone/unsigneddevicebuildpending.
REPORT-v42/README/VERIFICATION updatedhonestlyINPROGRESS. New291-workflowtimingtool
written+py_compileonly; fixedunfilteredfirst1000sourceprepared v42-1000-workflow-timing/.
DO NOTruntimebenchmarkuntilallownedCPUworkends; reference isRDKit215+FLINT2+
Python24+NumPy50, notRDKitkernel, no pinnedIPCnumericalcompatibilityclaim.

---

# CURRENT — host42 complete, Apple Simulator build active

All host42 jobs collected successfully: BulkIO58093 exit0; CSV70552 exit0;
CacheBulk68103 exit0; CacheIntegration17741 exit0.10000originalvalid, all291/
source/order/formula unchanged vs already C-ABI-verified fast41Swiftsnapshot,
CSV295columns2806165finitecells exact. Disk29112946bytes,reopen10000reuse0native.
Evidence v42-swift-10k/ andv42-host-integration-proof.json; logs archived.
One early CSVvalidator16755 exited1 missing not-yet-created catalog; diagnostic
retained, rerun after prerequisiteexit0 PASS. Not a chemistryfailure.
Host output tree safelyremoved afterall exits/proofs beforeAppletransition.
Full100k audit79422 stillactive using immutable42dylib outside removedhosttree;
latestparts0/1/3 at4160each zero differences, part2before2112. Do not alter
fivehashlocked reference tools/sourceCSV/dylib. Ownedprocesssamplepart0 showed
ordinary RDKit/Python reference BCUT work, notnativefactorrecursion.

Xcodegen72560 exit0collected. ACTIVE xcodebuild92355 Release build-for-testing
for owned iPad9AA7C21E..., log /private/tmp/chemob-v42-simulator-build.log.
OutputsONLY /private/tmp/chemob-descriptor-apple-target and
/private/tmp/chemob-descriptor-derived; nohosttree. Normal Japanese original
32D42BCF... remains booted with38; its6inputs backedup/hash/IDs beforeinstall,
app-v42/original-ipad-inputs-before.json andoriginal-ipad-before.json.
TestiPad9AA... andphoneACD... stillshutdown atbuildstart. Donotdelete/unboot
originalnormalSimulator. NewmetalUITestready;10kUITestaddsreview0/valid10000.
Needfocusedmetal+10kactualCSV+phone sample tests, originalPNGinspection,
normalJPinstallpreserve6, unsignedphysicalbuild, archive/cleanup byplatform.
Noactualphysicalspeed, noaskrepeat, noagents, norelease/commit/push.
REPORT-v42.md is explicitlyINPROGRESS; updateafterfullaudit/App gates.

---

# CURRENT — native42 all metal5 resolved; whole100k audit ACTIVE

Native42 frozen /private/tmp/chemob-descriptor-eval/libducksmiles-v42.dylib
SHA37243db55988cc222e2165069f0273f8cf106744b824cb143dee500018397ab2,
staticSHA19092fa06e83a6074ddf52aebac4aa5adb7998468d3fd144a5d359b9bea2c297.
Revisionducksmiles-workbench-v42-cluster291-rdkit2026036. Fullworkspace54680
ACTUALexit0:808passed/0failed/2ignored/0warnings, v42-workspace.json+greenlog.
All5 original metalclusters+permuted/synthetic158 all291+formulaPASS, zero
reference rejections. Native restricted canonical port matches423 raw rank-only
vectors including random graphs; TDD actualruntimeREDs/greenlogs retained.
Fullrank scope<=64/allbracket/allSingle/noStereo; NT one-donor prefix proof only
when uniquelygreatest initialclass. Never removesstereotags/inputorderguess.
Source-v42-native.tar.gz archive54339 exit0collected. Do not mutate native41/42
frozen dylibs or reference tools during active whole100k audit.

ACTIVE whole100k79422: parallel_catalog --information-backend flint --workers4,
output tasks/evidence/catalog289-sps/v42-full100k, log
/private/tmp/chemob-descriptor-v42-full100k.log. Four disjoint original ranges,
no native filtering. NOTYET whole100k completed/acceptanceclaim. Needcollectexit.
Host final42BulkIO58093 running original10k; cacheintegrationjuststarted,
cachebulk42Swiftbuild65216 pending. Check latesttool for sessionIDs.

Fixed41-factor2 actualSwift10k49368 exit0:10000valid, allnames/source/order,
CSVroundtrip, warmreuse10000/noadditionalnative. Independent295-column CSV
2806165finitecellsPASS. ExactoldDENSE41 vsfixedfactor2 vsactualSwift291 values
all10000PASS (validator92911 exit0); no tolerances. Evidence
v41-factor2-swift-10k/. Cold66.68s,warm123.74ms with otherownedjobs ->diagnostic
only, no fair speed/physicaliOSclaim. Largestisolatedpaired15.36x remainsvalid.

Real41negativeall5 Storecapture88773 exit0;fixture
apps/Chemob/Benchmarks/Fixtures/v41-clusters-unsupported-cache.json, provenance
v41-real-cluster-cache-provenance.json. FinalCacheIntegration expects42 plus
5arguments: old36cache,37cache,38cache,41clusternegative,remaining-metal5 JSON.
Addsall5recompute+5duplicates/diskreopen10reuse/all291exact.
NewiPad metalclusterUITestwritten but NOTexecuted. NormalinstalledAppstill38.
Needsfinal42host10k/cacheproof, archivehostwork,removeownedhostoutputtree then
AppleSim/UI/real10kCSV/phone/unsigneddevice/normalJP6update. Full100k audit uses
frozen libraryoutsidehosttree, so itsimmutableartifact cancontinue duringApple.

---

# CURRENT — guarded factor2 fully GREEN; metal-cluster references in progress

No native/App jobs active at checkpoint except reference generator27523
(tasks/descriptors/generate_cluster_reference.py, log /private/tmp/chemob-descriptor-cluster-reference.log).
All previous factor2 jobs collected: full Rust806passed/0failed/2ignored/0warnings,
source-v41-factor2-native.tar.gz archived, immutable factor2 dylibSHA
6b581740a27cd0e0730f0d7fad758afd7b0eb9776f177c448e9b3b160b433f14.
Actual paired measurement28611 exit0: median7.199659sec ->0.468792sec,
15.357896x; six WHOLE291-result JSON outputs byte-identical. MaxRSS
105840640 ->23183360bytes. MacM1 only, one598atom/632bond original molecule,
3alternating pairs fresh C processes, no cache/retry, no owned work overlap;
other host activity uncontrolled. v41-factor2-largest-timing/summary.json+run.log.
The older first-factor regression below remains REJECTED; never run oldBulkIO41
binary linked to bad variant. Rebuild against factor2 or final42.

Canonical non-stereo small-cluster Python partition prototype now passes all195
raw rank vectors (3original clusters x65atom permutations), evidence
cluster-canonical-prototype.json. Not production yet. Next generator creates
all5original/permuted all291 independent values plus random rank-only controls,
then actual RED before native port. Two NTclusters require only one donor with
unique highest initial canonical prefix, so no stereo-tag removal/rank guessing.
Current source41 prepared100k=99995/5unsupported. Full latest100k numerical
comparison, final App/10k/cache/UI/device build STILLpending. InstalledApp38.

---

# CRITICAL newest — first factor variant REJECTED; guarded version pending tests

Do NOTadoptfirst libducksmiles-v41-factor.dylib e82... despiteonecase16.8x.
Fullproductionworkspace73067 hung minutes in two oldstereo/SPS tests; native
sample showsdeep Solver::solve recursion. OwnedtestPID29843 terminatedSIGTERM;
Cargo73067 collected101, partialrunexplicitnotpass. RSS23MiB, no resourceleak;
unproductivecut expansion causedCPUregression. Evidence
v41-factor-unguarded-regression.json + terminatedlog/sample.
BulkIO41compiled69971collected0 but NEVEREXECUTED; it linksbadfirstfactorlibrary.
Do NOTrunthisbinary untilrebuiltwithnewfixedlibrary! No Appnewversioninstalled.

Actualshapeanalysis:85SPScontrolgraphs>32withcuts,59unbalanced; concrete128atom
macrocyclebestcutleaves120atoms, or120atomcageleaves118, causedmany repeated
largeproblem copies. cut-unbalanced-stereo-shapes.json retainssources/fixturelines.

AddedtworesourceguardsTDD. RuntimeRED62203 collected101,2expectedfailures:
ring48+48pendant graph shoulddeclineunproductivecut; explicitbudget1should
stoprecursivepath65. Thenpatchedfactorproduction: acceptcutONLY iflargest
component<=3n/4, countALL solvevisits(includingcachehits/unresolvedbranches),
defaultbudget4n+64; budgetNone propagates toexistingdensefallback. Exactmath/
revision41 unchanged. Currenttestjob inlatesttoolresult (cargo test --lib
--release information::), log /private/tmp/chemob-descriptor-cut-guard-green.log.
Collectactualexit beforeedits/build. No otheractivejobs. Needsnewlibrary
libducksmiles-v41-factor2.dylib (neveroverwritebadfirstfrozenlib), actualrevision
/hashassert, productionfullworkspace, largestpairedmeasurementRECHECK, latest
actual10kIO. Full100k numericauditNOTstarted yet (avoidbadfactorhang).

Correctoriginaldense41 baselinefull102+all291+formula PASS and800RusttestsPASS
remainvalid. Itsnative dlibSHA6b758b... immutable; normalAppstill38. Actual41
Storecachetestpassedoriginaldense41, notbadfirstfactorApp. Noagents/Applejobs.
Remainingfivecanonicalmetalcasesunimplemented. SourceprimaryRDKit hanoiSort.h
+FindRings.cpp downloaded withhashmanifest metal-canonical-next-sources.json;
newcanon.h/.cpp andMolOps.cpp alreadylocal. No metalrankportwrittenyet.

---

# LATEST — paired production speed proof PASS; fullworkspace job in latest tooloutput

Pairedmeasurement77887 ACTUALexit0 collected. Sixserializedresultfiles(all291)
BYTEIDENTICAL originaldense41 vsfactor41. Median11.509748sec ->0.683304sec =16.844x
on original598atom/632bond ChEBIsource. WholeCprocess maxRSS99729408 ->23166976bytes.
AppleM1 macOS,3alternatingpairs freshprocesses/no caches/nobufferretry. Noowned
compiler/auditoverlap; otherhostactivityuncontrolled; NOphysicaliOS/smallmolecule
generalclaim. summary+allrawoutputs v41-factor-largest-timing/. Userinformed.

Currentfull productionfactor workspace run starting NOW, logfile
/private/tmp/chemob-descriptor-v41-factor-workspace.log; sessionID latesttooloutput.
No sourceedits untilactualexit. Expected804pass/2ignored, NOTyetclaimed.
Next actualSwift10kIO viafrozenfactordylib (notmovingstaticwhilecargo runs),
retainedoriginaldataset, all10000nowexpectedvalidbutverify. ExistingBulkIO/main
usesChemistry.swift/Store + CoreSPMobjects. CheckexistingIOcommands/validator.
Needwhole100k numerical proof/finalAppand5metalcoordinationcasesstillpending.
No referenceaudits running currently. Native41-factordylibimmutable.

---

# Immediate — production factor41 built; paired measurement active

Frozen dense41 all102 original-large reference comparison25037 collectedexit0:
102accepted,all102all291+formula matching,zeroerrors. No tolerancechanges.
Proof v41-large102/comparison.json; logretained. CorecoverageupdatednewinputgatePASS
while whole100k latestnumericrerunstillpending; installednormalAppstill38.
Actual41 SwiftStore cachetest71903 collectedexit0, alloldcontrolsplusnewnegative
geometry/257C migration,reuse/reopen exactall291PASS; source41expectcompiledlinked
static ORIGINALdense41. Cachebuild/integrationlogsretained.

Nativefactor experimentalGREEN80117 collectedexit0:4tests,0.63sec, exactproducts
5000+i128/270digitsquare,58independent signedvectors,128randomblockgraphs against
independentdense recurrence, all599signedcoeffs of original598atom input.
AfterGREENproductionrouted n>32 viafactoring; strictadjacency/640/1280budget
validation extractedandapplied BEFORE traversal, <=32denseunchanged; biconnected
leaves densefallback, factorNone safelyfallbackdense. Extra invalid33vertex
/asymmetry/duplicate/641domaincontrolsadded. No runtime deps. Same mathematical
signedcoefficients/entropy reduction -> same semanticrevision41/cachevalidity;
newimplementationhash explicit, never overwritefrozenoriginaldense41.
Productionlibrarybuild82012 exit0collected, freshcclink+runtimeassert41 PASS:
/private/tmp/chemob-descriptor-eval/libducksmiles-v41-factor.dylib
SHAe82ecbb4cf21c9e8c2b041ed4a7a2be59a16906bf3347136a1a0556c9ba9995f
staticSHAe2de72d12fe270540fcbd83306ee362b7e23ccbcaad3665ce45e343b669eb003.
Metadata v41-factor-frozen-host-library.json. Fullproductionworkspace NOTyet run.

ACTIVE job77887 pairedlargestmeasurement only, logfile
/private/tmp/chemob-descriptor-v41-factor-largest-timing.log; source
tasks/descriptors/paired_large_measurement.py.3alternatingold/newpairs (6freshC
processes), original598atom/632bond SMILES, all291 + JSON, no cache/retry.
It requires WHOLEserializedoutputbytes identical all6 beforetimingclaim.
No ownedcompiler/reference/nativeaudit overlaps; otherhostactivityuncontrolled.
DO NOTstartCPUheavyjobs untilactualexit. Driver measure_one_native41 expects41;
olderCsource40retainedas evidence/v40-measure-one-native.c. No speedclaimyet.

AFTER pairedmeasurement collect0andinspectsummary, run fullRustworkspace and
whole100k numeric audit using FROZENfactor41 + explicitFLINT,4workers.
Reference toolscan'tchangewhileauditlocked. Compareeveryoriginal100kposition;
remain5ambiguousmetalunsupported,strictgatewillfail1expectedunlessresolvedfirst;
donotdeclarecomplete. Need stronger all102old/newbitparity orfinalwhole100krefs
beforeAppintegration; actualpairedlargestchecks onlyoneinput. All58+largest
coefficientproofsalreadyexact. Properchemical2D remainsfuture, no3D.

No Apple/Swift/Cargo jobsactive currently. Hosttargetstillpresent. Needfinal
engine/App actual10kIO/cache/UI/devicebuild andnormalJapaneseAppupdate preserving6.
Source-v41-native.tar.gz freezes originaldense41 BEFOREfactor; archivefinalvariant
afterfullGreen. Original38App/unsigneddeviceZIP/proofsremainseparate.

---

# Latest active — native factor experiment staged, NOT production

Full frozen41workspace16183 exit0 collected:800passed/2ignored/0warnings;
source-v41-native.tar.gz frozen BEFORE experimental factor module/stubs.
Cohort41parallel25037 STILLrunning, nofinalclaim; frozenlibraryimmutable,
reference toolshashlocked. Checklogbeforefinaluse.

New factor prototype58 exit0; largest598 prototype47240 exit0: all599signed
coefficients exactlymatch wholeFLINT;234memoizedsubgraphs,161leaves,maxleaf32,
73cutnodes. Originallargest independentlyretainedcoef fixture generation80576
exit0, information_largest_coefficients.tsv + metadata JSON.

Native experimentalmodule crates/smiles/src/information/factor.rs ONLY#[cfg(test)]
so productionstilloriginaldense41 algorithm. Integer::multiply ONLY#[cfg(test)].
Actual runtimeRED78527 collected101: two NotImplementedstub testsfail, asintended.
Then implemented exact signed multiplication, cut/componentpolynomial recursion,
balancedTarjan cuts, per-inputBTreeMapmemo, dense leaf fallback. Added128random
block graphs vs independentdense recurrence, plus599coefficient original598control.
Current testGREEN job80117 active: cargo test --lib --release factor_, log
/private/tmp/chemob-descriptor-cut-native-green.log. Collectactualexit/inspect.
No native productionrouting/performanceclaim yet! No sourceedits whilecompileractive.
Native coefficients() wrapper currently delegates to dense_coefficients(); factor
module uses dense_coefficients leaves to avoid future recursion. Need onlyafter
alltestsPASS wire moduleunconditionally and large(n>32)route AFTER original
strictadjacency/budgetvalidation. Preserveinvalidgraphs and n640/bonds1280guards.
New algorithm produces exact same signedcoeffvector and same entropyorder so
canretainsemanticrevision41 if provenbitidentical; newlibraryhash/variant explicit.
No Rustdependencies added. Need actualspeed/RSS proof after correctroute.

No other currentjobs beyond cohort25037 and factor80117. NoApple/Swift/process
benchmarksactive. Frozen38normalAppstillinstalled. Remainingfivecanonicalmetal
inputs stillunsupported, latest10kApp/whole100k numericalproof notyetdone.

---

# Immediate live state — native41 frozen; validation/prototype running

Actual41 build26722 collectedexit0, fresh cc linkexit0, Native runtime revision
assert41 PASS. Frozen /private/tmp/chemob-descriptor-eval/libducksmiles-v41.dylib
SHA6b758bee4d9af5815b8643952420c8f99001ed82f79e8a2972024972f3fbc9c8.
StaticSHA098f1a6e1a28eaf330dc45bdf9e5bafb9143ec1aa5bae4093913cf4c55833f8a.
Proof v41-frozen-host-library.json. Native source41 HARARY fix awaitingfulltests.

Active jobs:
* fullRustworkspace16183, /private/tmp/chemob-descriptor-v41-workspace.log.
* all102 original-large parallel25037 (fourworkers), log
/private/tmp/chemob-descriptor-v41-large102.log; evidence v41-large102/.
Lockedverificationtools cannotchange whilethisruns.
* NEW independent experiment47240 validates cut-polynomial decomposition on
maximum598original molecule against wholeFLINT. It may have poorperformance;
no nativeport yet and no speedclaim. Source prototype_cut_ipc.py, separate
from verificationoracle. All58 coefficient controls (forcedcutoff3) already
match, actual83357exit0collected, log chemob-descriptor-cut-ipc-prototype58.log.
Largestprototype finalstats shouldwrite cut-ipc-largest-prototype.json.
No need relaunch these jobs. Collect actualexit beforedependentactions.

No Apple/Swift/benchmarkjobs active. Cache41expectupdatedbutnotcompiled/run41.
InstallednormalApp38; do notrelabelfrozen38GUIproof. Full100k latestnumerical
andactual10knewengine/Appintegrationpending. Retain5unsupportedmetalcases.

---

# Current stage — Harary precision fix41 awaiting build/test, 2026-10-18

User resume active; do not finalize as complete. Native source now
`ducksmiles-workbench-v41-harary291-rdkit2026036`. Normal installed App STILL38.
No agents; Ryer catalog empty. No commit/push. Current host target
/private/tmp/chemob-descriptor-target (Swift scratch child swift-host).
Apple target/DerivedData/Generated platform trees absent. No physical phone tunnel.

## Live job (collect actual exit before dependencies/edits)
Cargo build native41 session26722; /private/tmp/chemob-descriptor-v41-library-build.log.
NO other active job. Then LINK FRESH dylib from STATIC archive via
cc -dynamiclib -Wl,-force_load,/private/tmp/chemob-descriptor-target/release/libducksmiles_smiles.a
-framework Security -framework CoreFoundation -liconv -o /private/tmp/chemob-descriptor-eval/libducksmiles-v41.dylib
Validate ACTUAL revision + hash using compare_catalog.Native before comparison.
Cargo manifest only staticlib+rlib: cargo build does NOT refresh old manually
linked release/libducksmiles_smiles.dylib! A v40 attempt copied stale39, got102
unsupported, excluded and retained INVALID under v40-stale-v39-library-rejected.
Correct40 linker/current revision9a2607... was subsequently validated.

## Frozen40 actual results
Prepared640atoms/1280bonds within4096bytes, raw sidecar256/512 unchanged.
Boundary/path257291value tests GREEN, fullRust799pass/2ignored/0warnings.
Whole100k prep99995accepted/5unsupported,102new/0regressions; evidencecomplete.
Frozen40 dylibSHA9a2607acd30540ab29ed41004e2cc118b03ab39ff28519d2e7dc82fc50bc23c5.
Native source snapshot source-v40-native.tar.gz frozen BEFORE new Harary test/fix.

Full102original-large cohort numerical comparison EXIT1 collected: all102accepted,
all formulas correct,95all291matching,7differ ONLY ds_Harary abs1.001e-9..1.439e-9,
unchangedtol1e-9. Currentcompleteproof v40-large102/comparison.json.
All102 pinned Ipc/AvgIpc numerical compatibility differ, SEPARATE from stable
integerdefinition agreement; no reference errors, no otherdescriptor mismatches.
Reference exactFraction distance-histogram sums for7 agree with denseNumPy values.
Fixture crates/smiles/tests/fixtures/large_harary.tsv; metadata
large-next-controls/harary-exact-reference.json. Genuine native runtimeRED
session63204 collected101, all7fail; evidence harary-runtime-red.log.
Production41 patch replaces naive reciprocal accumulation with exact distance
counts and fixed-order Kahan reduction, no tolerance change. Revision41 invalidates
oldresults. Need GREEN7, fullworkspace~800pass and new102parallel comparison.
Use --information-backend flint, --workers4 to avoid sequential7min reference
cohort bottleneck. FLINT exactoracle proven5tests + real132parallel CLI all291.
Do not change reference tools during active locked parallel comparison.

## Actual Swift cache work
CapturePreviousCache.swift utility compiled with actual frozen38 dylib, its
MachO original installname pointed to oldhost deps! install_name_tool changed
only captured EXECUTABLE dependency to frozen38; actual runtimeassert38passed.
Two actualStore unsupported records Pt@SP1 and257C captured and reopenedreuse2.
525bytefixture SHA901fc597648407e18e4c6d56b6d8c247610d049e1a1cc93d6219374d30fb2efb
Benchmarks/Fixtures/v38-geometry-large-unsupported-cache.json + provenanceevidence.
First invocation had extra --ui-testing and countassert trap133 BEFORE Store;
corrected utilityaccepts onlydirectory, explicitprivateTMPstore, green0.
CacheIntegration nowtakes THREE fixtureargs (argc4) v36overvalent/v37CX/v38new.
With40 actualStore testPASS: recomputed2, duplicate reuse2, disk reopenreuse4,
all291binary64 encode identical, no newcalls afterreload; priorcontrolsPASS.
Initial throwinprecondition compilefailure fixed by encoding outsideautoclosure.
Current source CacheIntegration revisionexpect41 (not yetcompiled/run41).
Host compiled40binary /private/tmp/chemob-descriptor-target/cache-integration40.
SPM Core5 objectfiles+Modules under swift-host/arm64-apple-macosx/release.

## Measured performance issue (needs concrete improvement)
Largest originalglobal52906/cohort51? DERIVE frompopulation102manifest,598atoms/632
bonds, originalinput1676bytes. C-only CLI with1MiB no outputretry all291 uncached
measured11–16sec, maxRSS~104MiB on loadedMac (Pythonreference/Swiftcompile overlap).
No physical iOS orfairpaired benchmarkclaim. Native first3sec sample2363/2363
samples in information::calculate/Integer::add_assign/memcpy/realloc. ExactIpc
coefficient dense recurrence is confirmed expensive, NOT simpleSMILES parsing.
C benchmark tasks/descriptors/measure_one_native.c currently expects40; retain
oldproof/source before using41. Fourcalls plus oneadditional/moresampled call
allfinished0, no active benchmarkjobs. Policy/files v40-largest-*.json/sample.
One sampler missed exitedprocess255, rerun with Python orchestration launched
sampler immediately, native0/sampler0 collected, actualprofile retained.
Python102 audit itself also spenttime in independent denseNumPy oracle, separate
from C engine timings. No claimsaboutfair speedratio.

After correctness41, consider mathematically exact articulation/component
polynomial factorization to avoid dense IPC on large sparse molecules, with
independent coefficient tests and exact old/new comparisons. No performance
implementation made yet. Current2fullmatrix store computeshalf+copiesotherhalf.
Integer already reuses Vecbuffers, no arbitrarylimb truncation.

Remaining5ambiguousmetal inputs requirepinnedcanonicalcoordination semantics:
18084,26896,53180,84673,98485. Diagnostic originals+reference ranks/datives retained.
Direct primary new_canon.h/.cpp + MolOps.cpp downloaded to
/private/tmp/chemob-descriptor-data/rdkit-{new_canon.h,new_canon.cpp,MolOps.cpp}.
Their hashes latesttoolhistory; canonical fullport NOT implemented. Native
parser canonical_ranks is NOT a compatible pinned rank substitute.

Need final whole100k numerical audit (frozen37 remainsbaseline sofar), latest
App/cache/actualdistinct10kIO/GUI/devicebuild integration. Expected first10k
input allsupported since remaining5globalindices>10000, but DO NOTclaim untilrun.
REPORT-v40.md scopesfailedcheckpoint. coverage latestinput checkpoint40failed7;
no goalcompletion. Goalbackend hadpaused despiteuserresume; do notcreateanother.
OriginalJapanese normal38Simulator preserves6savedIDs/name/source, leavebooted.

---

# Active stage — prepared640/1280 domain, v40 native pending population validation

Resume is active. No agents/commit/push. Ryer catalog empty. Current native
revision v40-large291-rdkit2026036; installed normal App remains v38.
Frozen39 adds valid @SP/@TB/@OH: whole100k preparation99893accepted/107unsupported,
zero regressions; original132geometry controls all291 match; full797pass/2ignored.
Evidence v39-acceptance100k.json and v39-geometry132/comparison.json.

v40 extends ONLY prepared workbench to640atoms/1280bonds; raw graph256/512 unchanged.
Actual runtimeRED two failed targets exit101; production patch then both GREEN0:
6 boundary controls include641atoms,1281bonds,8320bytes; path257 all291+formula
match independent reference, native runtime0.55sec (not fair speed benchmark).
No whole-population numerical40 claim yet. Whole100k PREPARATION job active:
session46850; log /private/tmp/chemob-descriptor-v40-acceptance100k.log.
It recompiles lib test due last comment-only source change. Collect actual exit
and success before cargo build/native freeze; do NOT mutate Rust inputs midbuild.
Host target /private/tmp/chemob-descriptor-target; Apple trees absent.

Reference-only python-flint0.9.0/FLINT3.6.0 installed pinned existing RDKit2026.03.6
env. Independent coefficient controls58+paths257/384/512/640 and validation/router
passed5tests24.609sec, actual exit0 collected. Exact Python oracle defaults to
old full-object recurrence <=256vertices/512edges and delegates larger to FLINT.
Optional --information-backend flint now transparently propagated by compare/
parallel_catalog.py, recorded in policy and module hashlocked. Native Rust
algorithm remains independent signed limbs/half matrices, no runtime FLINT dep.
Backend TDD red1 ->6tests green0; real parallel CLI132geometry,2workers all291+
formula green0. Logs/artifacts retained flint-cli-geometry132 and flint-* evidence.
No tolerance changed; pinned Ipc/AvgIpc numerical differences still separately
reported (stable field labels). Scripts changed since frozen37 tools archive.

All102 original >256atom unsupported cohort extracted independent of v40 status:
/private/tmp/chemob-descriptor-data/concrete/v40-large102.csv, SHA
1c6195712a5ea84d157eae19c749e05ea7d4c31ab5140d606da21949a11aa9bd.
Max598atoms/642bonds. Source positions+strings/provenance in
large-next-controls/population102-manifest.json. Need freeze new40 CABI after
successful cargo build, then all291 compare with explicitFLINT. Measure largest
actual native wall/RSS separately (Mac only) before performance redesign.

Remaining five small metal ambiguities indices18084,26896,53180,84673,98485;
remaining-metal5-diagnostic.json has original strings, reference ranks and dative
edges. Native parser canonical_ranks() is NOT pinned RDKit rank port; do NOT
substitute it arbitrarily. Later full rank policy/independent tests needed.
Current exact descriptors or strict100k goal NOT complete; no update_goal complete.
Backend goal state had paused despite user resume; no artificial new goal.
Original Japanese iPad normal38 installed; preserve original6saved IDs/name/source.
No Apple jobs active; owned test devices shutdown; phones unavailable DDI/tunnel.
Need eventual final engine/App/cache/10k/UI/devicebuild and full100k numeric proof.

---

# Active stage — non-tetrahedral TDD after frozen v38 App, 12:55 JST

User resume remains active; continue original accuracy/descriptors/App goal. No
agents. Ryer tool catalog rechecked EMPTY. No commit/push/release. Strict gate
still141unsupported, do NOT mark complete. Tsubame default, targeted rustfmt only.

## Current native source/test stage

Production native still v38 (NO geometry implementation yet). New integration
crates/smiles/tests/non_tetra.rs and two golden fixtures are added. Independent
ORIGINAL references now132positives/6invalid permutation inputs, all291+formula;
34population inputs plus synthetic permutation/charge/isotope/explicitH/arity
controls. Two source metal ambiguities18084/53180 explicitly excluded from this
focused cohort, retained in whole-source failed gate. Generator/data/metadata:
tasks/descriptors/generate_non_tetra_reference.py and evidence/catalog289-sps/
non-tetra-next-controls/. Metadata GetBoolProp(_ChiralityPossible) initiallyfailed
because the stored type is int; correctedGetIntProp, generator completesexit0.
FixtureSHA8cb448df1f42ce84a7380ba4223a22b33f74f8e7d644dd45340c6994a520163b.

Runtime RED127alreadycollected101:2testsPASS/positivegoldenFAIL allunsupported.
Latest132 runtimeRED RUNNING session45138? NO old45138 exited101collected. New
session ID is in latest tool result after this note; log /private/tmp/chemob-
descriptor-non-tetra132-runtime-red.log. COLLECT actualexit before productionedits.
Host target CURRENT /private/tmp/chemob-descriptor-target; Swift scratchchild
swift-host, no Apple trees. Current SwiftCore18host tests PASS, log+summaryarchived.
No active Swift/Apple/reference jobs. Only newnativeREDmayrun.

Planned implementation: prepared-only parse @SP/@TB/@OH into Atom.chirality, with
strict explicitpermutations1..3/1..20/1..30; omittednumber is valid. @SP0 etc
explicitzero INVALID. Negative/positive signs AFTER class are formalcharges, not
negativepermutations. PUBLIC raw parse/graph grammar must still reject these tags.
Add explicitallow flag to private parse_with_perception (threecalls); prepared
as-written true, rawnormalfalse. Newprivategeometryenum/helperonAtom is reasonable.

Pinneddirect source differs from stale webcache: actual Lipinski.cpp553linesSHA
427002f10cecccc044d8c0b4ab313ac467fccdae80e2d97488c6bda224fe51fd retained
/private/tmp/chemob-descriptor-data/rdkit-Lipinski-2026036.cpp. numAtomStereoCenters
counts presence of _ChiralityPossible; NumUnspecified ALSOrequires chiralTag==
CHI_UNSPECIFIED. Therefore do NOTexcludeNTfromalllegacypossiblecounts! A valid
NT tagwithtotaldegree4 may count1specified (CIP code can beNone). Addmetadata
flag forNumUnspecifiedwithoutinventingR/Scodes/refiningranks. ModernSPS onlycounts
Atom_Tetrahedral, so its eligible guard shouldskip retainedNTgeometries.

ConjugHybrid: validSP totaldegree2..4 overridesSP2D;TB2..5 ->SP3D;OH2..6 ->SP3D2.
Existing generic chiralitySome&&degree4 SP3 shortcut MUSTrestrict to @/@@. Add
SP2D enumvariant; charge.rs test-only exhaustivehybmatch needsSP2D label, runtime
charge matchalready has _ fallback. CleanupChirality clearsNTtag fordegreeoutside
class2..max, AFTER hybridization; it does NOTrejectthe molecule. New controls
include classoverspecarity and explicitHoverflow acceptedbyRDKit. Do NOTreject
these acceptedaritycases atparseguard. Before descriptors canclearinvalidNTtag;
itsfallbackhybridization is unchanged since outsideclassvaliddegree.

Restrict legacy implicit_hydrogen_transfers chirality predicate to regular @/@@,
otherwise incomplete neutralC@SP1H would incorrectly move bracketH toimplicit.
New132controls include this exactrisk. Current profile/SPS calculations operate
onoriginalpreparedgraph; no strippingannotations as productionnormalization.
Prototype/readtestall132+oldregressions beforecache39/Appintegration.

Direct pinned CPP sources (primarydownloaded;webcachewasstale) under
/private/tmp/chemob-descriptor-data/rdkit-{ConjugHybrid,Chirality,SmilesParse}-2026036.cpp;
SHAmanifest in next-controls/pinned-source-manifest.json. SeecleanupChirality in
Chirality.cpp. Legacy sourcehelperassign_code only@/@@; no falseCIPcodes forNT.

## Frozen v38 checkpoint DONE (installed normal App)

Native38 frozenSHA8cd3c5fcd76bab10c570f2f4b812a63640806b91c90ff0fe1201630057d9398a,
/private/tmp/chemob-descriptor-eval/libducksmiles-v38.dylib.217+24+50=291fields.
Original100k PREPARATION now99859accepted/141unsupported,107new/0regressions.
All107new ORIGINAL CX inputs CABI291+formula PASS;482CXgoldencontrolsPASS.
FullRust794pass/2ignored/0warnings. No whole100k numerical38rerunclaim: frozen37
whole99752accepted all291match baseline plus107incrementalproof separately.
Actual original10k Swift38 I/O:16unsupported, warm10000reuse/zeroNative;
291fields exactlymatch38CABI, previous37accepted all291 binary64unchanged.
Actual Simulator10kGUI import/search/export PASS. Its29132263byteCSV WHOLEFILE
byte-identicaltoMac;2801641finitecells and alloriginalsource/order/blankschecked.

VisualQA found roundednegativezero. iOSruntimeRED4assertions thenfixed helper;
4sharedformattertests+newCXUI allPASS. Core18host PASS. FinaliPhone17e flowPASS.
OriginalJapaneseiPad32D42BCF-E278-4EE5-AC74-D99C05D04963 updated normal38 (no testargs),
all6savedID/name/SMILES records preserved. Screenshot viewedoriginal; cLogP0.00.
Leave this originaldevice running! OwnedtestPad9AA7... andphoneACD3... shutdown0.
AllselectedoriginalPNGs/viewedmanifests/finalsuccess+REDlogs in app-v38/.

Finalunsignedphysicalarm64buildPASS,min17,families1/2,version0.1.0(1),fourCABI
symbols+three source noticesverified. ZIP addsrootMIT. Phonespairedbutunavailable
DDI/tunnel; NO physicalruntime/speedclaim. Deviceavailabilitysanitizedproofretained.
Frozen source-v38-final.tar.gz includesall5Rustworkspace/Coremanifest/App/project/
Coretests/public10kUITestfixtures+rights. Newer132NTtests are NOTin frozen38archive.
Report REPORT-v38.md; artifacts v38-checkpoint-artifacts.json. Archive38source is
frozen before nextnativeedits. Needlaterupdatefreshcurrentreports not relabeloldUI.
All ownedApplebuild/DerivedData/Generatedplatformtrees and rawxcresults/exports
removedAFTERactualexits+retainedselectedproofs. Receiptsapp-v38/. NoactiveApplejobs.

RemainingafterNT34support:102oversizedatoms(current256;max598,bondsmax642), five
ambiguousmetalpreparationinputs (3oldclusters+2NTclusters). Fullcanonicalranking
required bypinned metalcleanup, do NOTfakewithinputorder. Currentbestnextscope
isNT132controls first. Properchemical2D coordinates arefuture; no3D.

---

# Active delta — v38 native/host CX gate, 2026-10-18 12:20 JST

Do not confuse current source v38 with installed Simulator v37 or full numerical
100k frozen37 proof. Native v38 supports validated coordinate-free wU/wD ONLY,
retains source/cache identity, rejects unknown/radical/coordinate/wiggly blocks,
and excludes dative wedge edges. New module workbench_cx.rs. 482 original-reference
controls (107 population+375 synthetic) all291+formula PASS, runtime RED archived.
Full release workspace794pass/0fail/2ignored/0warnings. Whole100k preparation
now99859accepted/141unsupported,107new and ZERO regressions. New107 actualCABI
all291+formula against ORIGINAL annotated references PASS. Two pinned IPC numeric
differences among these107 retained separate. Do not claim current38 numerical
whole100k rerun; old37 all99752 numeric baseline +107 new verified separately.

Frozen38 /private/tmp/chemob-descriptor-eval/libducksmiles-v38.dylib SHA
8cd3c5fcd76bab10c570f2f4b812a63640806b91c90ff0fe1201630057d9398a.
Host tree CURRENT /private/tmp/chemob-descriptor-target, Swift scratchchild
swift-host. Core17tests allPASS after actual runtimeRED11assertions and whitespace
RED3; CompoundImport now preserves valid/unknown/malformed CX blocks before optional
SMI name; EngineInput permits spaces but blocks row injection. Two required
CacheIntegration args: v36-overvalent-cache.json AND v37-cx-unsupported-cache.json.
Actual frozen37 negative fixture captured with provenance;38 recomputes positive;
duplicate/diskrestore checks and all priorcache checks PASS. CacheIntegration,
BulkIO,CacheBulk executables compiled here. Actual original10k Swift38 BulkIO PASS:
16unsupported (nine newlyaccepted),291columns295CSVcols,source/order/reimport,warm
10000reuse/zero new native. JSON/CSV at /private/tmp/chemob-descriptor-data/concrete/
v38-bulk-io. Timing49s Mac under concurrent verification NOT physical/fair benchmark.
CSV validator just launched session36382; collect exit. Need all291 Swift/native
independent check via validate_swift_native_values.py. Apple source not yetbuilt38,
originalnormaluseriPad still37. Need final native features then Apple integration;
do not publish current38 UI proof from old37 screenshots.

Remaining non-tetra36 read-only experiment COMPLETED:34 candidate plain graphs
accepted, two ambiguousmetal clusters remain unsupported. Stripping tags changes
SPS on10 cases and NumUnspecifiedAtomStereoCenters on6. NEVER strip non-tetra
meaning as production support. Next consider prepared-only @SP/@TB/@OH retention,
correct hybridization class and modern/legacy stereo exclusion using pinnedrules.
No such production code has been written yet. Relevant pinned sources browsed:
ConjugHybrid.cpp setHybridization honors SP2Ddegree2..4/TB2..5/OH2..6;
MolOps.cpp organometallic cleanup uses FULL canonical ranks + changing dative counts
and atom processing order. Current coordination rejects >1metal candidate safely.
Do not fake ranking with input order. Three oldunclassified plus two Fe clusters
need this full policy; size102 remainsbound. Synthetic geometry controls/testRED
must precede code. Raw public parse/API should retain its prior grammar.

NO active nativebuild jobs; BulkIO process11874 exited0 collected. Prior Core81748
exit0 collected. CSVvalidator36382onlymaystillrun. Frozen107CABI comparer26119
exit0 collected. All checkpoint37 builds cleaned as documentedbelow. No agents.

---

# Active checkpoint — v37 frozen; unsupported-input investigation next

2026-10-18 12:10 JST. User resumed; no new stop. Keep the original objective:
RDKit-derived standard 2D definitions/accuracy + 50 additional nonduplicates,
App/API/cache/search/CSV, actual distinct 10k I/O, licensed source 100k validation.
No 3D. Do not mark complete while the required strict acceptance gate fails.
Goal backend still says paused despite user's resume; cannot set active via tools.
No budget, no new goal. No agents/commit/push/publish. Ryer catalog empty.

## Completed/current v37

Revision ducksmiles-workbench-v37-valence291-rdkit2026036; frozen dylib
/private/tmp/chemob-descriptor-eval/libducksmiles-v37.dylib
SHA 3e2fa996604baba34d8a8b0b38f5c77a64434fed18a2b4da7457cf16224cb882.
291 fields = 217 standards (215 pinned numerical, two stable exact IPC definitions)
+24 existing profiles +50 published extras. All source/API/App integrations done.

WHOLE 100k current v37 comparison COMPLETED: 99752 accepted, all 291 values and
formula within declared references/tolerances; ZERO numerical/reference/formula
errors. Overall strict acceptance gate FAILED for 248 valid unsupported source
rows. Pinned RDKit IPC numerical compatibility separately FAILED on 5617 rows.
Native exact IPC/AvgIPC independently passes all commonly accepted inputs;
relative/absolute policy explicit; CSV Ipc_stable/AvgIpc_stable prevents silent
model-feature substitution. Reference Python3.13.3/RDKit2026.03.6/NumPy2.5.4.
Comparison uses four original source ranges and no native acceptance filtering.
All processes exited/collected; earlier incomplete36 runs explicitly discarded.

Original 10k actual Swift + actual iPad Simulator I/O pass, 25 unsupported retained,
295 columns / 2799118 finite cells. Actual Simulator CSV all cells checked,
data-body byte-identical to Mac; headers use stable aliases. Latest v37 actual
ChemistryStore disk-capacity/reopen also PASS: 29042424 bytes /10000 restored /
zero new calculations, all values/source-order identical. Cache fixture uses real
v36 wrong-positive C5 and current37 rejects/recomputes. FullRust792pass/2ignored;
Core14pass;14603 synthetic valence acceptance/formula controls pass.

All latest iPad six-flow, actual10kexport, iPhone flow, focused actual visible
systemshare/cancel tests PASS. Three notices included, version0.1.0(1).
Physical generic arm64 iOS build PASS, unsigned ZIP retained. Phones unavailable
(DDI/tunnel unavailable), NO physical runtime/speed claims. Original Japanese
Simulator32D42BCF-E278-4EE5-AC74-D99C05D04963 left on latest normal v37 with all6
saved rows preserved. Fresh owned iPad9AA7C21E... and phoneACD3A432... shutdown,
exit0 collected. Never erase user original device or restore its old container.
Backbeat animation atlas/renderer already integrated, main Backbeat repo untouched.

## Current files/proofs

Report: tasks/evidence/catalog289-sps/REPORT.md (current, honest failed gates).
Full numeric proof: samefolder/v37-full100k/comparison.json and compressed diffs.
Cache: v37-cachebulk10k.json; App screenshots/logs/proofs: app-v37/.
Unsigned ZIP: app-v37/Chemob-0.1.0-unsigned-iphoneos.zip (MIT LICENSE alongside App).
Source: source-v37-final.tar.gz, includes all five workspace members, Core manifest,
App/resources/benchmarks and licenses; update archive docs if source changes.
Previous detailed operational history: RESUME-history-through-v37.md.

## Clean build state and rules

NO active build/reference/test jobs. Host target, Apple target, DerivedData,
Generated/iphoneos + iphonesimulator and owned raw xcresults/exports removed only
after actual exits, final summaries/selected original PNGs/logs retained. Receipts
are in current evidence. Source tree, installed original App, frozen external
libraries, source datasets/reference env retained.
One output tree per target; next native host slot is
/private/tmp/chemob-descriptor-target (currently absent). Swift scratch child
swift-host. Apple paths /private/tmp/chemob-descriptor-apple-target and
/private/tmp/chemob-descriptor-derived, remove/archive before transition.
Do not use cargo fmt --all. Targeted rustfmt skip_children=true only.
Tsubame indexed discovery default. No main duckSMILES/Backbeat repo edits.

## Next immediate work

Investigate the remaining 248 accepted-RDKit inputs without weakening correctness
or silently stripping stereochemical meaning. Read-only diagnostics retained:
v37-unsupported-diagnostics.json. 102 atom-limit (>256;8 also bond-limit>512),
108 CX annotations,36 non-tetrahedral tokens; flags overlap. Three unclassified
small metal clusters global26896/84673/98485. Largest excluded graph598atoms/642bonds.

First compare original annotated references against a frozen37 base-graph
experiment, WITHOUT changing production or claiming normalization equals support.
CX wedges can affect stereochemistry/SPS; never blindly strip extensions.
Use pinned Python /private/tmp/chemob-descriptor-env/bin/python (system lacks RDKit).
Original source/diffs preserve names/SMILES/indices. New implementation must start
with RED independent-reference tests and bump cache revision if results change.
Proper chemical 2D depiction is still future (current circle connectivity preview),
but not required for this descriptor milestone. Physical timing waits on device.

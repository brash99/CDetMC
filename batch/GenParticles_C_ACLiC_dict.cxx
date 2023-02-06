// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIhomedInandhudICDetOpticaldIbatchdIGenParticles_C_ACLiC_dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/home/nandhu/CDetOptical/batch/./GenParticles.C"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_GenParticles_C_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./GenParticles.C",
0
    };
    static const char* includePaths[] = {
"/usr/include/root",
"/usr/share/root",
"/usr/share/root/cling",
"/usr/include/root",
"/usr/include/root",
"/home/nandhu/CDetOptical/batch/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "GenParticles_C_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "GenParticles_C_ACLiC_dict dictionary payload"

#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "./GenParticles.C"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"GenParticles", payloadCode, "@",
"GenerateOneParticle", payloadCode, "@",
"InitOutput", payloadCode, "@",
"fE", payloadCode, "@",
"fIntRatio", payloadCode, "@",
"fM", payloadCode, "@",
"fMomFlatDist", payloadCode, "@",
"fMomMax", payloadCode, "@",
"fMomMean", payloadCode, "@",
"fMomMin", payloadCode, "@",
"fMomPowDist", payloadCode, "@",
"fOutFileName", payloadCode, "@",
"fP", payloadCode, "@",
"fPDG", payloadCode, "@",
"fPDGCode", payloadCode, "@",
"fPhiDist", payloadCode, "@",
"fPx", payloadCode, "@",
"fPy", payloadCode, "@",
"fPz", payloadCode, "@",
"fROOTFile", payloadCode, "@",
"fROOTTree", payloadCode, "@",
"fRand", payloadCode, "@",
"fThetaDist", payloadCode, "@",
"fVx", payloadCode, "@",
"fVy", payloadCode, "@",
"fVz", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("GenParticles_C_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_GenParticles_C_ACLiC_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_GenParticles_C_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_GenParticles_C_ACLiC_dict() {
  TriggerDictionaryInitialization_GenParticles_C_ACLiC_dict_Impl();
}

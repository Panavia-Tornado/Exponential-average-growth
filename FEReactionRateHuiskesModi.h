#pragma once
//=============================================================================
// This plugin example illustrates how to create a new material. 
// It requires FEBio 3.0 (or up)
//
// Author : Steve Maas
// Copyright (c) 2015 - 2020
// All rights reserved
//
//=============================================================================

//-----------------------------------------------------------------------------
// We need to include this file since our new material class will inherit from
// FEElasticMaterial, which is defined in this include file.

#pragma once
#include "FEBioMix/stdafx.h"
#include "FEBioMix/FEBiphasic.h"
#include "FEBioMix/FEMultiphasic.h"
#include <FEBioMech/FERemodelingElasticMaterial.h>
#include <FEBioMech/FEElasticMixture.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FECore/FEMesh.h>
#include "FEBioMix/FEChemicalReaction.h"
#include <FECore/FEMeshTopo.h>
#include <FECore/FEParameterList.h>

class FEReactionRateHuiskesModi : public FEReactionRate
{
public:
	//! constructor
	FEReactionRateHuiskesModi(FEModel* pfem);
	
    //! initialization
    virtual bool Init();
    
	//! reaction rate at material point
	virtual double ReactionRate(FEMaterialPoint& pt);
	
	//! tangent of reaction rate with strain at material point
	virtual mat3ds Tangent_ReactionRate_Strain(FEMaterialPoint& pt);
	
	//! tangent of reaction rate with effective fluid pressure at material point
	virtual double Tangent_ReactionRate_Pressure(FEMaterialPoint& pt);
	
public:
	FEParamDouble   m_B;					//!< mass supply coefficient
    FEParamDouble   m_psi0;					//!< specific strain energy at homeostasis
    double          m_D;                    //!< characteristic sensor distance

private:
    int             m_comp;                 //!< component of solid mixture (if applicable)
    std::vector<std::vector<int>>    m_EPL; //!< list of element proximity lists
    FEMeshTopo      m_topo;                 //!< mesh topology;
    bool            m_binit;                //!< initialization flag
    double          m_M;                    //!< molar mass of sbm
    int             m_lsbm;                 //!< local sbm value

	DECLARE_FECORE_CLASS();
};

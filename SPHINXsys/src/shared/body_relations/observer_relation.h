/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	observer_relation.h
 * @brief 	relation for observing.
 * @author	Xiangyu Hu
 */

#ifndef OBSERVER_RELATION_H
#define OBSERVER_RELATION_H

#include "contact_body_relation.h"
#include "observer_body.h"
namespace SPH
{
    /**
     * @class ObserverRelation
     * @brief The relation between an observer body and observed bodies
     */
    class ObserverRelation : public ContactRelation
    {
    public:
        ObserverRelation(ObserverBody &observer_body, RealBodyVector contact_bodies);
        virtual ~ObserverRelation(){};
    protected:
        bool is_configuration_updated_;
    };

}
#endif // OBSERVER_RELATION_H
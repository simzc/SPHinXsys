#include "observer_relation.h"

namespace SPH
{
    //=================================================================================================//
    ObserverRelation::ObserverRelation(ObserverBody &observer_body, RealBodyVector contact_bodies)
        : ContactRelation(observer_body, contact_bodies)
    {
        to_update_configuration_ = false;
        is_configuration_updated_ = false;
    }
    //=================================================================================================//
    void ObserverRelation::updateConfiguration()
    {
        if (!is_configuration_updated_)
        {
            ContactRelation::updateConfiguration();
        }
        is_configuration_updated_ = true;
    }
    //=================================================================================================//
    void ObserverRelation::updateRelation()
    {
        is_configuration_updated_ = false;
    }
    //=================================================================================================//
}

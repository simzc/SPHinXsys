#include "observer_relation.h"

namespace SPH
{
    //=================================================================================================//
    ObserverRelation(ObserverBody &observer_body, RealBodyVector contact_bodies)
        : ContactRelation(observer_body, contact_bodies),
          is_configuration_updated_(false)
    {
        to_update_configuration_ = false;
    }
    //=================================================================================================//
}

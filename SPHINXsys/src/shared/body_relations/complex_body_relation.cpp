#include "complex_body_relation.h"
#include "inner_body_relation.h"
#include "contact_body_relation.h"
#include "base_particle_dynamics.h"

namespace SPH
{
	//=================================================================================================//
	ComplexRelation::
		ComplexRelation(BaseInnerRelation &inner_relation, BaseContactRelation &contact_relation)
		: SPHRelation(inner_relation.sph_body_),
		  inner_relation_(inner_relation),
		  contact_relation_(contact_relation),
		  contact_bodies_(contact_relation_.contact_bodies_),
		  inner_configuration_(inner_relation_.inner_configuration_),
		  contact_configuration_(contact_relation_.contact_configuration_) {}
	//=================================================================================================//
	void ComplexRelation::updateConfigurationMemories()
	{
		inner_relation_.updateConfigurationMemories();
		contact_relation_.updateConfigurationMemories();
	}
	//=================================================================================================//
	void ComplexRelation::updateConfiguration()
	{
		inner_relation_.updateConfiguration();
		contact_relation_.updateConfiguration();
	}
	//=================================================================================================//
}

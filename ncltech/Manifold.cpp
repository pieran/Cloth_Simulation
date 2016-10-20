#include "Manifold.h"
#include "GameObject.h"
#include <nclgl\Matrix3.h>
#include <ncltech\NCLDebug.h>
#include "PhysicsEngine.h"

#define persistentThresholdSq 0.025f

typedef std::list<Contact> ContactList;
typedef ContactList::iterator ContactListItr;

Manifold::Manifold(PhysicsObject* nodeA, PhysicsObject* nodeB) : m_NodeA(nodeA), m_NodeB(nodeB)
{
}

Manifold::~Manifold()
{

}

void Manifold::ApplyImpulse()
{
	float softness = (m_NodeA->GetInverseMass() + m_NodeB->GetInverseMass()) / m_Contacts.size();
	for (Contact& contact : m_Contacts)
	{
		contact.normal.softness = softness;
		contact.normal.ApplyImpulse();

		//float friction_limit = sqrtf(2.0f) * contact.normal.impulseSum;
		float friction_limit = contact.normal.impulseSum;

		contact.friction1.impulseSumMin = -friction_limit;
		contact.friction1.impulseSumMax = friction_limit;
		
		contact.friction2.impulseSumMin = -friction_limit;
		contact.friction2.impulseSumMax = friction_limit;

		contact.friction1.ApplyImpulse();
		contact.friction2.ApplyImpulse();
	}
}

void Manifold::PreSolverStep(float dt)
{
	for (Contact& contact : m_Contacts)
	{
		UpdateConstraint(contact);
	}
}

void Manifold::UpdateConstraint(Contact& contact)
{

	Vector3 r1 = contact.relPosA;
	Vector3 r2 = contact.relPosB;


	Vector3 v1 = m_NodeA->GetLinearVelocity() + Vector3::Cross(m_NodeA->GetAngularVelocity(), r1);
	Vector3 v2 = m_NodeB->GetLinearVelocity() + Vector3::Cross(m_NodeB->GetAngularVelocity(), r2);


	Vector3 dv = v2 - v1;
	Vector3 tangent1 = dv - (contact.collisionNormal * Vector3::Dot(dv, contact.collisionNormal));
	if (Vector3::Dot(tangent1, tangent1) < 0.001f)
	{
		tangent1 = Vector3::Cross(contact.collisionNormal, Vector3(1, 0, 0));
		if (Vector3::Dot(tangent1, tangent1) < 0.001f)
		{
			tangent1 = Vector3::Cross(contact.collisionNormal, Vector3(0, 0, 1));
		}
	}
	tangent1.Normalise();

	Vector3 tangent2 = Vector3::Cross(contact.collisionNormal, tangent1);
	tangent2.Normalise();




	//Normal Collision Constraint
	float b = 0.0f;

	//Baumgarte Offset
	{
		const float baumgarte_scalar = 0.2f;
		const float baumgarte_slop = 0.01f;
		float penetration_slop = min(contact.collisionPenetration + baumgarte_slop, 0.0f);
		b += (baumgarte_scalar / PhysicsEngine::Instance()->GetDeltaTime()) * penetration_slop;
	}

	//Elasticity
	{
		const float elasticity = 0.8f;
		const float elasticity_slop = 0.5f;

		float elatisity_term = elasticity * Vector3::Dot(contact.collisionNormal,
			-m_NodeA->GetLinearVelocity()
			- Vector3::Cross(r1, m_NodeA->GetAngularVelocity())
			+ m_NodeB->GetLinearVelocity()
			+ Vector3::Cross(r2, m_NodeB->GetAngularVelocity())
			);

		b += min(elatisity_term + elasticity_slop, 0.0f);
	}


	//Friction Collision Constraints
	float friction = (m_NodeA->GetFriction() * m_NodeB->GetFriction());


	contact.normal = Constraint(m_NodeA, m_NodeB,
		-contact.collisionNormal,
		Vector3::Cross(-r1, contact.collisionNormal),
		contact.collisionNormal,
		Vector3::Cross(r2, contact.collisionNormal),
		b);

	contact.normal.impulseSumMin = 0.0f;

	contact.friction1 = Constraint(m_NodeA, m_NodeB,
		-tangent1 * friction,
		Vector3::Cross(-r1, tangent1) * friction,
		tangent1 * friction,
		Vector3::Cross(r2, tangent1) * friction,
		0.0f);

	contact.friction2 = Constraint(m_NodeA, m_NodeB,
		-tangent2 * friction,
		Vector3::Cross(-r1, tangent2) * friction,
		tangent2 * friction,
		Vector3::Cross(r2, tangent2) * friction,
		0.0f);

}

void Manifold::AddContact(const Vector3& globalOnA, const Vector3& globalOnB, const Vector3& normal, const float& penetration)
{
	Vector3 r1 = (globalOnA - m_NodeA->GetPosition());
	Vector3 r2 = (globalOnB - m_NodeB->GetPosition());


	Contact contact;
	contact.relPosA = r1;
	contact.relPosB = r2;
	contact.collisionNormal = normal;
	contact.collisionPenetration = penetration;

	m_Contacts.push_back(contact);
}

void Manifold::DebugDraw() const
{
	if (m_Contacts.size() > 0)
	{
		const Contact& c = m_Contacts.back();

		Vector3 globalOnA1 = m_NodeA->GetPosition() +  m_Contacts.back().relPosA;
		for (const Contact& contact : m_Contacts)
		{
			Vector3 globalOnA2 = m_NodeA->GetPosition() + contact.relPosA;
			Vector3 globalOnB = m_NodeB->GetPosition() + contact.relPosB;

			NCLDebug::DrawThickLine(globalOnA1, globalOnA2, 0.02f, Vector4(0.0f, 1.0f, 0.0f, 1.0f));
			NCLDebug::DrawPoint(globalOnA2, 0.05f, Vector4(0.0f, 0.5f, 0.0f, 1.0f));

			NCLDebug::DrawThickLine(globalOnB, globalOnA2, 0.01f, Vector4(1.0f, 0.0f, 1.0f, 1.0f));

			globalOnA1 = globalOnA2;
		}
	}
}
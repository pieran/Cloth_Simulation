

#pragma once

#include "NCLDebug.h"
#include "PhysicsEngine.h"

class Spring
{
public:
	Spring(PhysicsObject* objA, PhysicsObject* objB,
		const Vector3& globalOnA, const Vector3& globalOnB)
	{
		this->objA = objA;
		this->objB = objB;

		Vector3 ab = globalOnB - globalOnA;
		this->distance = ab.Length();
		k_stiffness = 0.001f;
		k_damping = 0.01f;

		Vector3 r1 = (globalOnA - objA->GetPosition());
		Vector3 r2 = (globalOnB - objB->GetPosition());

		this->localOnA = Matrix3::Transpose(objA->GetOrientation().ToMatrix3()) * r1;
		this->localOnB = Matrix3::Transpose(objB->GetOrientation().ToMatrix3()) * r2;
	}

	void ApplyForce()
	{
		Vector3 globalOnA = objA->GetWorldSpaceTransform() * localOnA;
		Vector3 globalOnB = objB->GetWorldSpaceTransform() * localOnB;

		Vector3 ab = globalOnB - globalOnA;
		Vector3 abn = ab;
		abn.Normalise();



		Vector3 r1 = objA->GetOrientation().ToMatrix3() * localOnA;
		Vector3 r2 = objB->GetOrientation().ToMatrix3() * localOnB;

		Vector3 globalOnA = r1 + objA->GetPosition();
		Vector3 globalOnB = r2 + objB->GetPosition();

		Vector3 ab = globalOnB - globalOnA;
		Vector3 abn = ab;
		abn.Normalise();

		this->j1 = -abn;
		this->j2 = Vector3::Cross(-r1, abn);
		this->j3 = abn;
		this->j4 = Vector3::Cross(r2, abn);

		//Baumgarte Offset (Adds energy to the system to counter slight solving errors that accumulate over time - known as 'constraint drift')
		{
			float distance_offset = distance - ab.Length();
			float baumgarte_scalar = 0.1f;
			b = -(baumgarte_scalar / PhysicsEngine::Instance()->GetDeltaTime()) * distance_offset;
		}


		delta = 0.0f;

		//J * M(-1) * J(t)
		float contstraint_mass = objA->GetInverseMass() * Vector3::Dot(j1, j1)
			+ Vector3::Dot(j2, (objA->GetInverseInertia() * j2))
			+ objB->GetInverseMass() * Vector3::Dot(j3, j3)
			+ Vector3::Dot(j4, (objB->GetInverseInertia() * j4))
			+ softness;

		if (contstraint_mass > 0.00001f)
		{
			//JV
			float jv = Vector3::Dot(j1, objA->GetLinearVelocity())
				+ Vector3::Dot(j2, objA->GetAngularVelocity())
				+ Vector3::Dot(j3, objB->GetLinearVelocity())
				+ Vector3::Dot(j4, objB->GetAngularVelocity());

			float denom = -(jv + b);
			delta = denom / contstraint_mass;
		}



		float oldImpulseSum = impulseSum;
		impulseSum = min(max(impulseSum + delta, impulseSumMin), impulseSumMax);
		float realDelta = impulseSum - oldImpulseSum;

		objA->SetLinearVelocity(objA->GetLinearVelocity() + (j1 * realDelta) * objA->GetInverseMass());
		objA->SetAngularVelocity(objA->GetAngularVelocity() + objA->GetInverseInertia() * (j2 * realDelta));
		objB->SetLinearVelocity(objB->GetLinearVelocity() + (j3 * realDelta) * objB->GetInverseMass());
		objB->SetAngularVelocity(objB->GetAngularVelocity() + objB->GetInverseInertia() * (j4 * realDelta));
	}

	virtual void DebugDraw() const
	{
		Vector3 globalOnA = objA->GetOrientation().ToMatrix3() * localOnA + objA->GetPosition();
		Vector3 globalOnB = objB->GetOrientation().ToMatrix3() * localOnB + objB->GetPosition();

		NCLDebug::DrawThickLine(globalOnA, globalOnB, 0.02f, Vector4(0.0f, 0.0f, 0.0f, 1.0f));
		NCLDebug::DrawPoint(globalOnA, 0.05f, Vector4(1.0f, 0.8f, 1.0f, 1.0f));
		NCLDebug::DrawPoint(globalOnB, 0.05f, Vector4(1.0f, 0.8f, 1.0f, 1.0f));
	}

protected:
	PhysicsObject* objA;
	PhysicsObject* objB;

	float   distance;
	float	k_stiffness;
	float	k_damping;

	Vector3 localOnA, localOnB;
};
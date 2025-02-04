#include "KSTrajAdiabaticSpinDerivative.h"

namespace Kassiopeia
{

KSTrajAdiabaticSpinDerivative::KSTrajAdiabaticSpinDerivative() = default;
KSTrajAdiabaticSpinDerivative::KSTrajAdiabaticSpinDerivative(const KSTrajAdiabaticSpinDerivative&) = default;
KSTrajAdiabaticSpinDerivative::~KSTrajAdiabaticSpinDerivative() = default;

void KSTrajAdiabaticSpinDerivative::AddToTime(const double& aTime)
{
    fData[0] += aTime;
    return;
}
void KSTrajAdiabaticSpinDerivative::AddToSpeed(const double& aSpeed)
{
    fData[1] += aSpeed;
    return;
}
void KSTrajAdiabaticSpinDerivative::AddToVelocity(const KGeoBag::KThreeVector& aVelocity)
{
    fData[2] += aVelocity.X();
    fData[3] += aVelocity.Y();
    fData[4] += aVelocity.Z();
    return;
}
void KSTrajAdiabaticSpinDerivative::AddToForce(const KGeoBag::KThreeVector& aForce)
{
    fData[5] += aForce.X();
    fData[6] += aForce.Y();
    fData[7] += aForce.Z();
    return;
}
void KSTrajAdiabaticSpinDerivative::AddToMDot(const double& anMDot)
{
    fData[8] += anMDot;
}
void KSTrajAdiabaticSpinDerivative::AddToPhiDot(const double& aPhiDot)
{
    fData[9] += aPhiDot;
}
}  // namespace Kassiopeia

/*
 * KStaticElectroMagnetField.hh
 *
 *  Created on: 25 Mar 2016
 *      Author: wolfgang
 */

#ifndef KSTATICELECTROMAGNETFIELD_HH_
#define KSTATICELECTROMAGNETFIELD_HH_

#include "KMagneticFieldSolver.hh"
#include "KMagnetostaticField.hh"
#include "KSmartPointer.hh"

namespace KEMField
{

class KStaticElectromagnetField : public KMagnetostaticField
{
  public:
    KStaticElectromagnetField();
    ~KStaticElectromagnetField() override;

    void SetDirectory(const std::string& aDirectory);
    void SetFile(const std::string& aFile);

    void SetFieldSolver(KSmartPointer<KMagneticFieldSolver> solver);
    KSmartPointer<KMagneticFieldSolver> GetFieldSolver();

    void SetContainer(KSmartPointer<KElectromagnetContainer> aContainer);
    KSmartPointer<KElectromagnetContainer> GetContainer() const;

  protected:
    void InitializeCore() override;
    void CheckSolverExistance() const;

    KThreeVector MagneticPotentialCore(const KPosition& aSamplePoint) const override;
    KThreeVector MagneticFieldCore(const KPosition& aSamplePoint) const override;
    KGradient MagneticGradientCore(const KPosition& aSamplePoint) const override;
    std::pair<KThreeVector, KGradient> MagneticFieldAndGradientCore(const KPosition& P) const override;

  private:
    KSmartPointer<KElectromagnetContainer> fContainer;
    KSmartPointer<KMagneticFieldSolver> fFieldSolver;

    std::string fFile;
    std::string fDirectory;
};

} /* namespace KEMField */

#endif /* KSTATICELECTROMAGNETFIELD_HH_ */

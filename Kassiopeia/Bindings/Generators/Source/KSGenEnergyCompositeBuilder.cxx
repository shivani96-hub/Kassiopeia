#include "KSGenEnergyCompositeBuilder.h"
#include "KSGenValueFixBuilder.h"
#include "KSGenValueSetBuilder.h"
#include "KSGenValueListBuilder.h"
#include "KSGenValueUniformBuilder.h"
#include "KSGenValueBoltzmannBuilder.h"
#include "KSGenValueGaussBuilder.h"
#include "KSGenValueGeneralizedGaussBuilder.h"
#include "KSRootBuilder.h"

#ifdef Kassiopeia_USE_ROOT
#include "KSGenValueFormulaBuilder.h"
#include "KSGenValueHistogramBuilder.h"
#endif

using namespace Kassiopeia;
using namespace std;

namespace katrin
{

    template< >
    KSGenEnergyCompositeBuilder::~KComplexElement()
    {
    }

    STATICINT sKSGenEnergyCompositeStructure =
        KSGenEnergyCompositeBuilder::Attribute< string >( "name" ) +
        KSGenEnergyCompositeBuilder::Attribute< string >( "energy" ) +
        KSGenEnergyCompositeBuilder::ComplexElement< KSGenValueFix >( "energy_fix" ) +
        KSGenEnergyCompositeBuilder::ComplexElement< KSGenValueSet >( "energy_set" ) +
        KSGenEnergyCompositeBuilder::ComplexElement< KSGenValueList >( "energy_list" ) +
        KSGenEnergyCompositeBuilder::ComplexElement< KSGenValueUniform >( "energy_uniform" ) +
        KSGenEnergyCompositeBuilder::ComplexElement< KSGenValueBoltzmann >( "energy_boltzmann" ) +
        KSGenEnergyCompositeBuilder::ComplexElement< KSGenValueGauss >( "energy_gauss" ) +
        KSGenEnergyCompositeBuilder::ComplexElement< KSGenValueGeneralizedGauss >( "energy_generalized_gauss" ) ;

    STATICINT sKSGenEnergyComposite =
        KSRootBuilder::ComplexElement< KSGenEnergyComposite >( "ksgen_energy_composite" );

#ifdef Kassiopeia_USE_ROOT
    STATICINT sKSGenEnergyCompositeStructureROOT =
            KSGenEnergyCompositeBuilder::ComplexElement< KSGenValueFormula >( "energy_formula" ) +
            KSGenEnergyCompositeBuilder::ComplexElement< KSGenValueHistogram >( "energy_histogram" );
#endif

}

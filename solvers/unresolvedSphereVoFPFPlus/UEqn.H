MRF.correctBoundaryVelocity(U);

volScalarField rhoAlpha("rhoAlpha", rho*coupling.alpha());

fvVectorMatrix UEqn
(
    fvm::ddt(rhoAlpha, U)
  - fvm::Sp(fvc::ddt(rhoAlpha), U)
  + fvm::div(rhoPhi, U)
  - fvm::Sp(fvc::div(rhoPhi), U)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(rho, U)
  + fvm::Sp(coupling.Sp(), U)
 ==
    - coupling.Su()
  + fvOptions(rho, U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        coupling.alpha()*
        fvc::reconstruct
        (
            (
                mixture.surfaceTensionForce()
              - ghf*fvc::snGrad(rho)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    fvOptions.correct(U);
}

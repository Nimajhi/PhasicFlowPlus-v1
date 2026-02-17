{
    if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volScalarField alpha("alpha", coupling.alpha());
    surfaceScalarField alphaf("alphaf", fvc::interpolate(alpha));

    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA*alpha)
      + MRF.zeroFilter(fvc::interpolate(alpha*rho*rAU())*fvc::ddtCorr(U, phi, Uf))
    );
    MRF.makeRelative(phiHbyA);

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*alphaf*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    surfaceScalarField rAUfAlpha("rAUfAlpha", rAUf*alphaf);
    surfaceScalarField rAUfAlpha2("rAUfAlpha2", rAUf*alphaf*alphaf);

    constrainPressure(p_rgh, U, phiHbyA, rAUfAlpha, MRF);

    while (pimple.correctNonOrthogonal())
    {
        volScalarField ddtVoidFraction = fvc::ddt(alpha);

        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUfAlpha2, p_rgh) == fvc::div(phiHbyA) + ddtVoidFraction
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(p_rgh.select(pimple.finalInnerIter()));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUfAlpha);
            U.correctBoundaryConditions();

            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    fvc::correctUf(Uf, U, phi);
    fvc::makeRelative(phi, U);

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    if (!correctPhi)
    {
        rAU.clear();
    }
}

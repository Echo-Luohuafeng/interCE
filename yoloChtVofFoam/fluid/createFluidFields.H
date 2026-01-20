PtrList<twoPhaseMixtureEThermo> thermoFluid(fluidRegions.size());

PtrList<temperaturePhaseChangeTwoPhaseMixture> mixtureFluid(fluidRegions.size());

// 基本场
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<surfaceScalarField> rhoPhiFluid(fluidRegions.size()); // VOF 需要的质量通量

// 压力与重力相关
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());

// 湍流模型: VOF 通常使用不可压缩湍流基类
PtrList<incompressible::turbulenceModel> turbulenceFluid(fluidRegions.size());

// VOF 特有场
PtrList<interfaceProperties> interfaceFluid(fluidRegions.size()); // 界面属性
PtrList<volScalarField> kappaEffFluid(fluidRegions.size());       // 有效热导率
PtrList<volScalarField> rhoCpFluid(fluidRegions.size());          // 体积热容



// MRF 和 fvOptions
PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

// 压力参考点 (用于 p_rgh)
List<label> pRefCellFluid(fluidRegions.size());
List<scalar> pRefValueFluid(fluidRegions.size());

// 重力场 (全局唯一)
const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);


// 2. 遍历所有流体区域进行初始化
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    // --- 初始化 p_rgh (先于 U 以满足某些依赖，或保持一致性) ---
    Info<< "    Adding to p_rghFluid\n" << endl;
    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    // --- 初始化 U ---
    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    // --- 初始化 phi ---
    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(UFluid[i]) & fluidRegions[i].Sf()
        )
    );

    // --- 初始化 VOF 热力学模型 (twoPhaseMixtureEThermo) ---
    Info<< "    Adding to thermoFluid (twoPhaseMixtureEThermo)\n" << endl;
    thermoFluid.set
    (
        i,
        new twoPhaseMixtureEThermo(UFluid[i], phiFluid[i])
    );

    // --- 初始化相变模型 (Mixture) ---
    Info<< "    Adding to mixtureFluid\n" << endl;
    mixtureFluid.set
    (
        i,
        temperaturePhaseChangeTwoPhaseMixture::New(thermoFluid[i], fluidRegions[i])
    );

    // --- 初始化密度 rho (基于相分数混合) ---
    Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            // rho = alpha1*rho1 + alpha2*rho2
            thermoFluid[i].alpha1()*thermoFluid[i].rho1() 
          + thermoFluid[i].alpha2()*thermoFluid[i].rho2()
        )
    );
    rhoFluid[i].oldTime(); // 求解 ddt(rho) 需要旧时间步值

    // --- 初始化质量通量 rhoPhi ---
    Info<< "    Adding to rhoPhiFluid\n" << endl;
    rhoPhiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "rhoPhi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(rhoFluid[i])*phiFluid[i]
        )
    );

    // --- 初始化界面属性 (Interface) ---
    Info<< "    Adding to interfaceFluid\n" << endl;
    interfaceFluid.set
    (
        i,
        new interfaceProperties
        (
            thermoFluid[i].alpha1(),
            UFluid[i],
            thermoFluid[i].transportPropertiesDict()
        )
    );

    // --- 初始化湍流模型 (Incompressible) ---
    Info<< "    Adding to turbulenceFluid\n" << endl;
    turbulenceFluid.set
    (
        i,
        incompressible::turbulenceModel::New
        (
            UFluid[i],
            phiFluid[i],
            thermoFluid[i]
        ).ptr()
    );

    // --- 处理重力项 gh 和 ghRef (用于 p_rgh 计算) ---
    Info<< "    Adding to hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(word::null, dimLength, Zero)
        )
    );

    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    ghFluid.set(i, new volScalarField("gh", (g & fluidRegions[i].C()) - ghRef));
    ghfFluid.set(i, new surfaceScalarField("ghf", (g & fluidRegions[i].Cf()) - ghRef));


    // --- 强制 p_rgh 与 p 一致 ---
    // CHT 标准操作，但在 VOF 中 p 通常是主变量。
    // 这里我们参照 createFields.H 中的逻辑：p_rgh = p - rho*gh
    volScalarField& p = thermoFluid[i].p();

    // 读取 PIMPLE 字典用于参考压力设置
    const dictionary& pimpleDict = fluidRegions[i].solutionDict().subDict("PIMPLE");
    
    pRefCellFluid[i] = -1;
    pRefValueFluid[i] = 0.0;

    if (p_rghFluid[i].needReference())
    {
        setRefCell
        (
            p,
            p_rghFluid[i],
            pimpleDict,
            pRefCellFluid[i],
            pRefValueFluid[i]
        );
    }
    
    // 如果 p_rgh 是读取的，更新 p；反之亦然。
    // 在这里我们先更新 p_rgh 以匹配 p 和 rho (注意 rho 刚被计算)
    p_rghFluid[i] = p - rhoFluid[i]*ghFluid[i];
    
    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());
    fluidRegions[i].setFluxRequired(thermoFluid[i].alpha1().name());


    // --- 初始化 kappaEff (有效热导率) ---
    Info<< "    Adding to kappaEffFluid\n" << endl;
    kappaEffFluid.set
    (
        i,
        new volScalarField
        (
            IOobject("kappaEff", runTime.timeName(), fluidRegions[i], IOobject::NO_READ, IOobject::NO_WRITE),
            thermoFluid[i].kappa() // 初始值为层流 kappa，循环中会加入湍流部分
        )
    );

    // --- 初始化 rhoCp (体积热容) ---
    Info<< "    Adding to rhoCpFluid\n" << endl;
    rhoCpFluid.set
    (
        i,
        new volScalarField
        (
            IOobject("rhoCp", runTime.timeName(), fluidRegions[i], IOobject::NO_READ, IOobject::NO_WRITE),
            rhoFluid[i]*thermoFluid[i].Cp()
        )
    );
    rhoCpFluid[i].oldTime();



    // --- MRF 和 fvOptions ---
    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set(i, new IOMRFZoneList(fluidRegions[i]));

    Info<< "    Adding fvOptions\n" << endl;
    fluidFvOptions.set(i, new fv::options(fluidRegions[i]));

    turbulenceFluid[i].validate();
}

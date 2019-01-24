within ;
package TankSystem

  model DigiLab_TankLevel_Flat

    extends Modelica.Icons.Example;
    import g = Modelica.Constants.g_n;

    // Tank related variables and parameters
    parameter Real qInflow(unit="m3/s") = 0.1     "Flow through input valve";
    parameter Real area(unit="m2") = 1            "Cross sectional area of the tank";
    parameter Real areaOut(unit="m2") = 0.04      "Cross sectional area of outflow tube";
    parameter Real h0(unit="m") = 0               "Initial filling level";
    Real           h(start=h0,unit="m")           "Tank filling level";
    Real           qOutflow(unit="m3/s")          "Flow through output valve";

    // Controller related variables and parameters
    parameter Real K_p = 0.1                "p-Gain";
    parameter Real K_i = 1                  "i-Gain";
    parameter Real K_d = 1                  "d-Gain";
    parameter Real T(unit="s") = 1          "Time constant of the tank system";
    parameter Real h_setpoint = 1           "Setpoint level for control";
    Real           opening(min=0, max=1)      "Valve opening in the range [0-1]";

    Modelica.Blocks.Continuous.LimPID PID(
      k = K_p,
      Ti = K_i,
      Td = K_d,
      yMax = 1,
      yMin = 0)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));


  equation

    assert(h>=0, "Tank filling level must be greater or equal to zero");

    der(h) = (qInflow - qOutflow)/area;   // Mass balance equation
    qOutflow = opening * areaOut * sqrt(2 * g * h);

    PID.u_s = h;
    PID.u_m = h_setpoint;
    opening = PID.y;

  annotation (
      experiment(
        StopTime=180,
        __Dymola_NumberOfIntervals=1000,
        Tolerance=1e-06));
  end DigiLab_TankLevel_Flat;

  package Obsolete
    extends Modelica.Icons.ObsoleteModel;

    model FlatTank

      extends Modelica.Icons.Example;

      // Tank related variables and parameters
      parameter Real flowLevel(unit="m3/s")=0.02;
      parameter Real area(unit="m2") = 1            "cross sectional area of the Tank";
      parameter Real flowGain(unit="m2/s")= 0.05;
      Real           h(start=0,unit="m")            "Tank level";
      Real           qInflow(unit="m3/s")           "Flow through input valve";
      Real           qOutflow(unit="m3/s")          "Flow through output valve";

      // Controller related variables and parameters
      parameter Real K_p=2                          "p-Gain";
      parameter Real K_i=2                          "i-Gain";
      parameter Real K_d=2                          "d-Gain";
      parameter Real T(unit="s")= 10                "Time constant of the tank system";
      parameter Real minV=0.0001, maxV=10;          // Limits for flow output
      parameter Real ref = 0.25                     "Reference level for control";
      Real           error                          "Deviation from reference level";
      Real           outCtr                         "Control signal without limiter";
      Real           x                              "State variable x for controller";
      Real           y                              "State variable y for controller";

      // Limiter function to reflect min. and max. flows through output valve
      encapsulated function limitValue
        input Real pMin;
        input Real pMax;
        input Real p;
        output Real pLim;
      algorithm
        pLim := if p>pMax then pMax
                else if p<pMin then pMin
                else p;
      end limitValue;

    equation

      assert( minV>=0,"minV must be greater or equal to zero");

      der(h) = (qInflow-qOutflow)/area;   // Mass balance equation
      qInflow  = flowLevel;
      qOutflow = FlatTank.limitValue(
            minV,
            maxV,
            -flowGain*outCtr);
      error  = ref-h;
      der(x) = error/T;
      y = T*der(error);
      outCtr = K_p*error+K_i*x+K_d*y;

    annotation (
        experiment(StartTime = 0, StopTime = 250, Tolerance = 1e-6, Interval = 0.5));
    end FlatTank;

    model DigiLab_TankLevel_Flat_error

      extends Modelica.Icons.Example;
      import g = Modelica.Constants.g_n;

      // Tank related variables and parameters
      parameter Real qInflow(unit="m3/s") = 0.2    "Flow through input valve";
      parameter Real area(unit="m2") = 1            "cross sectional area of the tank";
      parameter Real areaOut(unit="m2") = 0.05      "cross sectional area of outflow tube";
      Real           h(start=0,unit="m")            "Tank filling level";
      Real           qOutflow(unit="m3/s")          "Flow through output valve";

      // Controller related variables and parameters
      parameter Real K_p = 2                        "p-Gain";
      parameter Real K_i = 2                        "i-Gain";
      parameter Real K_d = 2                        "d-Gain";
      parameter Real T(unit="s") = 10               "Time constant of the tank system";
      parameter Real setpoint = 0.25                "Setpoint level for control";
      Real           error                          "Deviation from reference level";
      Real           outCtr(start=0)                "Control signal without limiter";
      Real           x(start=0)                     "State variable x for controller";
      Real           y                              "State variable y for controller";
      Real           alpha;

    equation

      assert(h>=0, "Tank filling level must be greater or equal to zero");

      der(h) = (qInflow - qOutflow)/area;   // Mass balance equation
      qOutflow = alpha * areaOut * sqrt(2 * g * h);
      error  = (-1) * (setpoint - h);

      if outCtr >= 1 or outCtr <= 0 then
        der(x) = 0;
      else
        der(x) = error/T;
      end if;

      y = T*der(error);
      outCtr = max(0, min(1, K_p * (error + K_i*x + K_d*y)));
      alpha = outCtr;

    annotation (
        experiment(StartTime = 0, StopTime = 250, Tolerance = 1e-6, Interval = 0.5));
    end DigiLab_TankLevel_Flat_error;

    model Tank

      extends Modelica.Icons.Example;

      inner Modelica.Fluid.System system(allowFlowReversal=false)
        annotation (Placement(transformation(extent={{70,50},{90,70}})));
      Modelica.Fluid.Valves.ValveIncompressible valve(
        redeclare package Medium =
            Modelica.Media.Water.ConstantPropertyLiquidWater,
        m_flow_nominal=20,
        dp_nominal=1000)
        annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
      Modelica.Fluid.Vessels.OpenTank tank(
        nPorts=2,
        redeclare package Medium =
            Modelica.Media.Water.ConstantPropertyLiquidWater,
        crossArea=1,
        level_start=0,
        use_portsData=false,
        height=2)
        annotation (Placement(transformation(extent={{-60,0},{-20,40}})));
      Modelica.Fluid.Sources.MassFlowSource_h inlet(
        nPorts=1,
        redeclare package Medium =
            Modelica.Media.Water.ConstantPropertyLiquidWater,
        m_flow=20)
        annotation (Placement(transformation(extent={{-90,-30},{-70,-10}})));
      Modelica.Fluid.Sources.Boundary_pT outlet(nPorts=1, redeclare package
          Medium = Modelica.Media.Water.ConstantPropertyLiquidWater)
        annotation (Placement(transformation(extent={{90,-30},{70,-10}})));
      Modelica.Blocks.Continuous.LimPID PID(
        yMax=1,
        yMin=0,
        controllerType=parametrs.controllerType,
        k=parametrs.k,
        Ti=parametrs.Ti,
        Td=parametrs.Td)
               annotation (Placement(transformation(extent={{40,10},{20,30}})));
      Modelica.Blocks.Sources.RealExpression setpoint(y=1)
        annotation (Placement(transformation(extent={{60,-10},{40,10}})));
      Modelica.Blocks.Sources.RealExpression sensor(y=tank.level)
        annotation (Placement(transformation(extent={{-30,30},{-10,50}})));
      Parametrs parametrs(
        controllerType=Modelica.Blocks.Types.SimpleController.PID,
        k=1,
        Ti=10,
        Td=1) annotation (Placement(transformation(extent={{-90,50},{-70,70}})));
    equation
      connect(inlet.ports[1], tank.ports[1]) annotation (Line(points={{-70,-20},{
              -44,-20},{-44,0}}, color={0,127,255}));
      connect(valve.port_a, tank.ports[2]) annotation (Line(points={{-10,-20},{
              -36,-20},{-36,0}}, color={0,127,255}));
      connect(valve.port_b, outlet.ports[1])
        annotation (Line(points={{10,-20},{70,-20}}, color={0,127,255}));
      connect(PID.y, valve.opening)
        annotation (Line(points={{19,20},{0,20},{0,-12}}, color={0,0,127}));
      connect(sensor.y, PID.u_s) annotation (Line(points={{-9,40},{50,40},{50,20},
              {42,20}}, color={0,0,127}));
      connect(setpoint.y, PID.u_m)
        annotation (Line(points={{39,0},{30,0},{30,8}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=250));
    end Tank;

    record Parametrs
      extends Modelica.Icons.Record;
      import Modelica.SIunits;
      import Modelica.Blocks.Types.SimpleController;

      parameter SimpleController controllerType=SimpleController.PID "Type of controller";
      parameter Real k(min=0, unit="1") = 1 "Gain of controller";
      parameter SIunits.Time Ti(min=Modelica.Constants.small)=0.5 "Time constant of Integrator block"
       annotation(Dialog(enable=controllerType==SimpleController.PI or controllerType==SimpleController.PID));
      parameter SIunits.Time Td(min=0)= 0.1 "Time constant of Derivative block"
       annotation(Dialog(enable=controllerType==SimpleController.PD or controllerType==SimpleController.PID));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Parametrs;

    model Tank_Image

      extends Modelica.Icons.Example;

      Modelica.Fluid.Valves.ValveIncompressible Valve(
        redeclare package Medium =
            Modelica.Media.Water.ConstantPropertyLiquidWater,
        m_flow_nominal=20,
        dp_nominal=1000)
        annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
      OpenTank_Image Tank(
        nPorts=1,
        redeclare package Medium =
            Modelica.Media.Water.ConstantPropertyLiquidWater,
        crossArea=1,
        level_start=0,
        use_portsData=false,
        height=2)
        annotation (Placement(transformation(extent={{-60,0},{-20,40}})));
      Modelica.Fluid.Sources.MassFlowSource_h Inlet(
        redeclare package Medium =
            Modelica.Media.Water.ConstantPropertyLiquidWater,
        m_flow=20,
        nPorts=1)
        annotation (Placement(transformation(extent={{-90,50},{-70,70}})));
      Modelica.Fluid.Sources.Boundary_pT Outlet(nPorts=1, redeclare package
          Medium = Modelica.Media.Water.ConstantPropertyLiquidWater)
        annotation (Placement(transformation(extent={{90,-30},{70,-10}})));
      Modelica.Blocks.Continuous.LimPID PID(
        controllerType=Modelica.Blocks.Types.SimpleController.PID,
        yMax=1,
        yMin=0)
               annotation (Placement(transformation(extent={{40,10},{20,30}})));
      Modelica.Blocks.Sources.RealExpression Setpoint(y=1)
        annotation (Placement(transformation(extent={{60,-10},{40,10}})));
      Modelica.Blocks.Sources.RealExpression Sensor(y=Level)
        annotation (Placement(transformation(extent={{-30,30},{-10,50}})));
    equation
      connect(Valve.port_a,Tank. ports[1]) annotation (Line(points={{-10,-20},{
              -40,-20},{-40,0}}, color={0,127,255}));
      connect(Valve.port_b,Outlet. ports[1])
        annotation (Line(points={{10,-20},{70,-20}}, color={0,127,255}));
      connect(PID.y,Valve. opening)
        annotation (Line(points={{19,20},{0,20},{0,-12}}, color={0,0,127}));
      connect(Sensor.y, PID.u_s) annotation (Line(points={{-9,40},{50,40},{50,20},
              {42,20}}, color={0,0,127}));
      connect(Setpoint.y, PID.u_m)
        annotation (Line(points={{39,0},{30,0},{30,8}}, color={0,0,127}));
      connect(Inlet.ports[1], Tank.port_a) annotation (Line(points={{-70,60},{-60,
              60},{-60,40}}, color={0,127,255}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=250));
    end Tank_Image;

  model OpenTank_Image "Simple tank with inlet/outlet ports"
      import Modelica.Constants.pi;

    // Tank properties
    Modelica.SIunits.Height level(stateSelect=StateSelect.prefer, start=
            level_start_eps) "Level height of tank";
    Modelica.SIunits.Volume V(stateSelect=StateSelect.never) "Actual tank volume";

    // Tank geometry
    parameter Modelica.SIunits.Height height "Height of tank";
    parameter Modelica.SIunits.Area crossArea "Area of tank";

    // Ambient
    parameter Medium.AbsolutePressure p_ambient=system.p_ambient
        "Tank surface pressure"
      annotation(Dialog(tab = "Assumptions", group = "Ambient"));
    parameter Medium.Temperature T_ambient=system.T_ambient
        "Tank surface Temperature"
      annotation(Dialog(tab = "Assumptions", group = "Ambient"));

    // Initialization
    parameter Modelica.SIunits.Height level_start(min=0)=0.5*height
        "Start value of tank level" annotation (Dialog(tab="Initialization"));

    // Mass and energy balance, ports
    extends Modelica.Fluid.Vessels.BaseClasses.PartialLumpedVessel(
      final fluidVolume = V,
      final fluidLevel = level,
      final fluidLevel_max = height,
      final vesselArea = crossArea,
      heatTransfer(surfaceAreas={crossArea+2*sqrt(crossArea*pi)*level}),
      final initialize_p = false,
      final p_start = p_ambient);

      Modelica.Fluid.Interfaces.FluidPort_a port_a annotation (Placement(
            transformation(extent={{-120,100},{-100,120}}), iconTransformation(
              extent={{-106,94},{-94,106}})));
    protected
    final parameter Modelica.SIunits.Height level_start_eps=max(level_start,
          Modelica.Constants.eps);

  equation
    // Total quantities
    V = crossArea*level "Volume of fluid";
    medium.p = p_ambient;

    // Source termsEnergy balance
      if Medium.singleState or energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState then
      Wb_flow = 0
          "Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)";
    else
      Wb_flow = -p_ambient*der(V);
    end if;

    //Determine port properties
    for i in 1:nPorts loop
      vessel_ps_static[i] = max(0, level - portsData_height[i])*system.g*medium.d + p_ambient;
    end for;

  initial equation
      if massDynamics == Modelica.Fluid.Types.Dynamics.FixedInitial then
      level = level_start_eps;
      elseif massDynamics == Modelica.Fluid.Types.Dynamics.SteadyStateInitial then
      der(level) = 0;
    end if;

      annotation (defaultComponentName="tank",
        Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.2), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={255,255,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.VerticalCylinder),
            Line(points={{-100,100},{-100,-100},{100,-100},{100,100}}, color={0,0,
                  0})}),
        Documentation(info="<html>
<p>
Model of a tank that is open to the ambient at the fixed pressure
<code>p_ambient</code>.
</p>
<p>
The vector of connectors <b>ports</b> represents fluid ports at configurable heights, relative to the bottom of tank.
Fluid can flow either out of or in to each port.
</p>
The following assumptions are made:
<ul>
<li>The tank is filled with a single or multiple-substance medium having a density higher than the density of the ambient medium.</li>
<li>The fluid has uniform density, temperature and mass fractions</li>
<li>No liquid is leaving the tank through the open top; the simulation breaks with an assertion if the liquid level growths over the height.</li>
</ul>
<p>
The port pressures represent the pressures just after the outlet (or just before the inlet) in the attached pipe.
The hydraulic resistances <code>portsData.zeta_in</code> and <code>portsData.zeta_out</code> determine the dissipative pressure drop between tank and port depending on
the direction of mass flow. See <a href=\"modelica://Modelica.Fluid.Vessels.BaseClasses.VesselPortsData\">VesselPortsData</a> and <i>[Idelchik, Handbook of Hydraulic Resistance, 2004]</i>.
</p>
<p>
With the setting <code>use_portsData=false</code>, the port pressure represents the static head
at the height of the respective port.
The relationship between pressure drop and mass flow rate at the port must then be provided by connected components;
Heights of ports as well as kinetic and potential energy of fluid entering or leaving are not taken into account anymore.
</p>
</html>",   revisions="<html>
<ul>
<li><i>Dec. 12, 2008</i> by Ruediger Franke: move port definitions
   to BaseClasses.PartialLumpedVessel; also use energy and mass balance from common base class</li>
<li><i>Dec. 8, 2008</i> by Michael Wetter (LBNL):<br>
Implemented trace substances.</li>
<li><i>Jan. 6, 2006</i> by Katja Poschlad, Manuel Remelhe (AST Uni Dortmund),
   Martin Otter (DLR):<br>
   Implementation based on former tank model.</li>
<li><i>Oct. 29, 2007</i> by Carsten Heinrich (ILK Dresden):<br>
Adapted to the new fluid library interfaces:
<ul> <li>FluidPorts_b is used instead of FluidPort_b (due to it is defined as an array of ports)</li>
    <li>Port name changed from port to ports</li></ul>Updated documentation.</li>
<li><i>Apr. 25, 2006</i> by Katrin Pr&ouml;l&szlig; (TUHH):<br>
Limitation to bottom ports only, added inlet and outlet loss factors.</li>
</ul>
</html>"));
  end OpenTank_Image;
  end Obsolete;
  annotation (uses(Modelica(version="3.2.2")));
end TankSystem;

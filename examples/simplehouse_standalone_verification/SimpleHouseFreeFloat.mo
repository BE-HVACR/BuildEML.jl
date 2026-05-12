within ;
model SimpleHouseFreeFloat

  package MediumAir = Buildings.Media.Air;


  parameter Modelica.Units.SI.Area AWall=100 "Wall area";
  parameter Modelica.Units.SI.Area AWin=5 "Window area";
  parameter Real gWin(min=0, max=1, unit="1") = 0.3 "Solar heat gain coefficient of window";
  parameter Modelica.Units.SI.Volume VZone=AWall*3 "Wall area";
  parameter Modelica.Units.SI.HeatFlowRate QHea_flow_nominal=700
    "Nominal capacity of heating system";
  parameter Modelica.Units.SI.MassFlowRate mWat_flow_nominal=QHea_flow_nominal/
      10/4200 "Nominal mass flow rate for water loop";
  parameter Modelica.Units.SI.MassFlowRate mAir_flow_nominal=VZone*2*1.2/3600
    "Nominal mass flow rate for air loop";

  parameter Modelica.Units.SI.PressureDifference dpAir_nominal=200
    "Pressure drop at nominal mass flow rate for air loop";
  parameter Boolean allowFlowReversal=false
    "= false because flow will not reverse in these circuits";

  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor walCap(T(fixed=true),
      C=10*AWall*0.05*1000*1000)
    "Thermal mass of walls"
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=270,
        origin={158,-16})));
  Buildings.Fluid.MixingVolumes.MixingVolume
                                   zon(
    redeclare package Medium = MediumAir,
    V=VZone,
    nPorts=0,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    m_flow_nominal=mAir_flow_nominal,
    massDynamics=Modelica.Fluid.Types.Dynamics.DynamicFreeInitial)
    "Very based zone air model"
    annotation (Placement(transformation(extent={{148,34},{168,14}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor conRes(R=1/2/AWall)
    "Thermal resistance for convective heat transfer with h=2" annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={98,4})));
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3
                                            weaDat(filNam=
       "C:/Users/rzhsong/Downloads/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.mos")
    "Weather data reader"
    annotation (Placement(transformation(extent={{-192,-26},{-172,-6}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor walRes(R=0.25/AWall/0.04)
              "Thermal resistor for wall: 25 cm of rockwool"
    annotation (Placement(transformation(extent={{48,-26},{68,-6}})));
  Buildings.HeatTransfer.Sources.PrescribedTemperature TOut
    "Exterior temperature boundary condition"
    annotation (Placement(transformation(extent={{-92,-26},{-72,-6}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow win
    "Very simple window model"
    annotation (Placement(transformation(extent={{48,-66},{68,-46}})));
  Modelica.Blocks.Math.Gain gaiWin(k=AWin*gWin)
    "Gain for window solar transmittance and area as HGloHor is in W/m2"
    annotation (Placement(transformation(extent={{8,-66},{28,-46}})));
  Buildings.BoundaryConditions.WeatherData.Bus
                                     weaBus "Weather data bus"
    annotation (Placement(transformation(extent={{-152,-26},{-132,-6}}),
        iconTransformation(extent={{-160,-10},{-140,10}})));
equation
  connect(conRes.port_a,zon. heatPort)
    annotation (Line(points={{98,14},{98,24},{148,24}},     color={191,0,0}));
  connect(weaDat.weaBus,weaBus)  annotation (Line(
      points={{-172,-16},{-142,-16}},
      color={255,204,51},
      thickness=0.5));
  connect(walRes.port_b,walCap. port) annotation (Line(points={{68,-16},{148,-16}},
                         color={191,0,0}));
  connect(TOut.T,weaBus. TDryBul)
    annotation (Line(points={{-94,-16},{-118,-16},{-118,-15.95},{-141.95,-15.95}},
                                                          color={0,0,127}));
  connect(TOut.port,walRes. port_a)
    annotation (Line(points={{-72,-16},{48,-16}}, color={191,0,0}));
  connect(gaiWin.y,win. Q_flow) annotation (Line(points={{29,-56},{48,-56}},
                           color={0,0,127}));
  connect(gaiWin.u,weaBus. HGloHor) annotation (Line(points={{6,-56},{-141.95,-56},
          {-141.95,-15.95}},     color={0,0,127}));
  connect(conRes.port_b,walCap. port) annotation (Line(points={{98,-6},{98,-16},
          {148,-16}},                                   color={191,0,0}));
  connect(win.port,walCap. port) annotation (Line(points={{68,-56},{98,-56},{98,
          -16},{148,-16}},                                 color={191,0,0}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-211.75,44},{-32.25,-76}},
          fillColor={238,238,238},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None),
        Rectangle(
          extent={{-12,44},{188,-76}},
          fillColor={238,238,238},
          fillPattern=FillPattern.Solid,
          pattern=LinePattern.None),
        Text(
          extent={{52.5,24.5},{-16.5,43.5}},
          textColor={0,0,127},
          fillColor={255,213,170},
          fillPattern=FillPattern.Solid,
          textString="Building"),
        Text(
          extent={{-153,25},{-211,43}},
          textColor={0,0,127},
          fillColor={255,213,170},
          fillPattern=FillPattern.Solid,
          textString="Weather")}),
    uses(Buildings(version="11.0.0"), Modelica(version="4.0.0")));
end SimpleHouseFreeFloat;

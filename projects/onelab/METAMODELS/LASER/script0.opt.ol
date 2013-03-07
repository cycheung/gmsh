//MAX TEMP AT DERMIS-EPIDERMIS INTERFACE
Plugin(Probe).View = 0;
Plugin(Probe).X=1.e-6;
Plugin(Probe).Y=OL.eval(OL.get(Parameters/Skin/3DERMIS)*1.e-3);
Plugin(Probe).Z=0;
Plugin(Probe).Run;
Plugin(MinMax).View=1;
Plugin(MinMax).OverTime=1;
Plugin(MinMax).Argument=0;
Plugin(MinMax).Run;
Save View [3] "tempMaxInterface.txt"; 
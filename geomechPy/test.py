import Utilities as UT



UT = UT.Utilities()

#[dtco,dtsm] = UT.slowness2velocity(200,300,"km/s")
pressure=UT.convert_pressure(1000,'GPa','Pa')

print(pressure)
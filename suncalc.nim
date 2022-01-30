import pylib
import times, math
const rad = PI / 180
# ported from suncalc-py.
#[ discard LICENSE: """
suncalc-py is ported from suncalc.js under the BSD-2-Clause license.
Copyright (c) 2014, Vladimir Agafonkin
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:
  1. Redistributions of source code must retain the above copyright notice, this list of
    conditions and the following disclaimer.
  2. Redistributions in binary form must reproduce the above copyright notice, this list
    of conditions and the following disclaimer in the documentation and/or other materials
    provided with the distribution.
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  """
  ]#
type SunCoords = object
  dec*: float64
  ra*: float64
type Position* = object
  azimuth*: float64
  altitude*: float64
# sun times configuration (angle, morning name, evening name)
type SunTime = object
  angle*:float64
  morningName*, eveningName*:string
  morningTime*:DateTime
  eveningTime*:DateTime

proc initSunTime*(n: float64, before, after: string): SunTime =
  result.angle = n
  result.morningName = before
  result.eveningName = after
let DEFAULT_TIMES: seq[SunTime] = @[
  initSunTime(-0.833, "sunrise", "sunset"),
  initSunTime(-0.3, "sunrise_end", "sunset_start"),
  initSunTime(-6, "dawn", "dusk"),
  initSunTime(-12, "nautical_dawn", "nautical_dusk"),
  initSunTime(-18, "night_end", "night"),
  initSunTime(6, "golden_hour_end", "golden_hour")
]
# date/time constants and conversions
const dayMs = 1000 * 60 * 60 * 24
const J1970 = 2440588
const J2000 = 2451545

proc to_julian(date: DateTime): auto =
  date.toTime().toUnix * 1000 / dayMs - 0.5 + J1970

proc from_julian(j: float64): DateTime =
  let ms_date = int (j + 0.5 - J1970) * dayMs
  let dur = initDuration(milliSeconds = ms_date)
  return initDateTime(1, mJan, 1970, 0, 0, 0) + dur

proc to_days(date: DateTime): float64 =
  to_julian(date) - J2000
# general calculations for position
# obliquity of the Earth
let e = rad * 23.4397

proc right_ascension(l, b: float64): float64 = arctan2((sin(l) * cos(e) - tan(b) *
    sin(e)), cos(l))

proc declination(l, b: float64): float64 = arcsin(sin(b) * cos(e) + cos(b) *
    sin(e) * sin(l))

proc azimuth(H, phi, dec: float64): float64 = arctan2(sin(H), cos(H) * sin(phi) -
    tan(dec) * cos(phi))

proc altitude(H, phi, dec: float64): float64 =
  return arcsin(sin(phi) * sin(dec) + cos(phi) * cos(dec) * cos(H))

proc sidereal_time(d, lw: float64): float64 =
  return rad * (280.16 + 360.9856235 * d) - lw

proc solar_mean_anomaly(d: auto): auto =
  return rad * (357.5291 + 0.98560028 * d)

proc ecliptic_longitude(M: auto): auto =
  # equation of center
  let C = rad * (1.9148 * sin(M) + 0.02 * sin(2 * M) + 0.0003 * sin(3 * M))
  # perihelion of the Earth
  let P = rad * 102.9372
  return M + C + P + PI

proc sun_coords(d: auto): SunCoords =
  let M = solar_mean_anomaly(d)
  let L = ecliptic_longitude(M)
  return SunCoords(dec: declination(L, 0), ra: right_ascension(L, 0));
# calculations for sun times
const J0 = 0.0009

proc julian_cycle(d, lw: auto): auto =
  return round(d - J0 - lw / (2 * PI))

proc approx_transit(Ht, lw, n: float64): float64 =
  # echo "approx:",Ht," ", lw," ", n," ", J0," ", PI
  return J0 + (Ht + lw) / (2 * PI) + n

proc solar_transit_j(ds, M, L: auto): auto =
  return J2000 + ds + 0.0053 * sin(M) - 0.0069 * sin(2 * L)

proc hour_angle(h, phi, d: float64): float64 =
  result = arccos((sin(h) - sin(phi) * sin(d)) / (cos(phi) * cos(d)))
  assert not result.isNaN

proc observer_angle(height: int): auto =
  return -2.076 * sqrt(height.float) / 60

proc get_set_j(h:float64, lw, phi, dec, n, M, L: float64): auto =
  discard """Get set time for the given sun altitude
  """
  let w = hour_angle(h, phi, dec)
  assert not w.isNaN
  let a = approx_transit(w, lw, n)
  return solar_transit_j(a, M, L)

proc get_position*(date:DateTime, lng, lat: auto): Position =
  discard """Calculate sun position for a given date and latitude/longitude
  """
  let
    lw = rad * -lng
    phi = rad * lat
    d = to_days(date)
    c = sun_coords(d)
    H = sidereal_time(d, lw) - c.ra
  return Position(azimuth: azimuth(H, phi, c.dec),
    altitude: altitude(H, phi, c.dec))
# import times

proc addOffset*(dt:DateTime, offset:Duration = initDuration(seconds = -dt.utcOffset())) : DateTime =
  # echo "offset from utc:" & $(-dt.utcOffset())
  result = dt + offset

proc get_time*(date:DateTime, lng, lat: auto, height:int = 0,
    sunAngle: SunTime): SunTime = 
  let
    lw = rad * -lng
    phi = rad * lat
    dh = observer_angle(height)
    d = to_days(date)
    n = julian_cycle(d, lw)
    ds = approx_transit(0, lw, n)
    M = solar_mean_anomaly(ds)
    L = ecliptic_longitude(M)
    dec = declination(L, 0)
    Jnoon = solar_transit_j(ds, M, L)
  let
    h0 : float64 = (sunAngle.angle + dh) * rad
    Jset = get_set_j(h0, lw, phi, dec, n, M, L)
    # Jrise = Jnoon - (Jset - Jnoon)
    Jrise = 2 * Jnoon - Jset
  let jrise = from_julian(Jrise)
  let jset = from_julian(Jset)
  var s: SunTime = sunAngle
  s.morningTime = jrise# .addOffset()# + Duration(seconds: jrise.utcOffset())
  s.eveningTime = jset# .addOffset()
  return s
let sunRiseAngle = initSunTime(-0.833, "sunrise", "sunset")

proc sunRise*(date:DateTime, lng, lat: auto, height:int = 0) :DateTime = 
  result = get_time(date, lng,lat,height, sunRiseAngle).morningTime

proc sunSet*(date:DateTime, lng, lat: auto, height:int = 0) :DateTime= 
  result = get_time(date, lng,lat,height, sunRiseAngle).eveningTime

proc get_times*(date:DateTime, lng, lat: auto, height:int = 0,
    sunAngles: seq[SunTime] = DEFAULT_TIMES): seq[SunTime] =
  for time in sunAngles:
    # date:DateTime, lng, lat: auto, height:int = 0, sunAngle: SunTime
    result.add get_time(date,lng,lat,height,time)
  return result

# workaround: daylight savings already changed.
let baseOffset*: Duration = initDuration(seconds = -1 * initDateTime(1,mJan,2000,1,0,0,0).utcOffset())

proc example() =
  ## 
  let date = now()
  let lng = 34.801
  let lat = 31.867
  for i in 0..100:
    var rise = sunRise(date + initDuration(days = i)  , lng, lat)
    # can use rise.utcOffset() which returns # of seconds, but times are wrong by 1 hour 
    # (daylight savings time was accounted for twice)
    # so I use a base offset from Jan 1st 
    var truerise = rise.addOffset(baseOffset)
    echo truerise


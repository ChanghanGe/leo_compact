#
# * PROPOSED LARGE BROADBAND CONSTELLATIONS - 2010s MEGACONSTELLATIONS
# *
# * Starlink
# *
# * This is the 'perfect' fully-deployed first-generation SpaceX
# * Starlink system. Only the initial "first shell" of the planned
# * initial Ku/Ka-band constellations is simulated here; this is
# * likely to be the "inner shell" for quite some time. Higher
# * shells and the lower, later, V-band constellations are not yet
# * simulated.
# *
# * Initial single launch of 60 satellites in May 2019, after
# * two test satellites were launched in February 2018. Further
# * launches have been carried out.
# *
# * This was based on the November 2018 FCC filing, which modified
# * and slightly reduced the first shell authorised in March 2018,
# * as well as lowering the satellites to 550km from 1,150km altitude.
# * The planned geometry of further added shells, if any, can also be
# * expected to be altered.
# *
# * The planned design of this first shell changed again in August
# * 2019, when SpaceX advised the FCC of a redistribution of the
# * orbital planes and satellites, while still keeping the same
# * overall number of satellites.
# *
# * Attempting to simulate any further shells does not appear
# * worthwhile, especially since SpaceX filed for an additional
# * 30,000 satellites in October 2019. That would have led to
# * a multiple-shell constellation of over 42,000 satellites.
# * Though SaVi could simulate that many satellites, little
# * would be learned from doing so. In April 2019 SpaceX filed
# * again with the FCC, lowering the altitude of all proposed
# * shells and satellites to under 575km, but leaving this first
# * shell unchanged from the earlier August 2019 filing.
# *
# * Expected user terminal minimum elevation angle for this shell
# * is 25 degrees, rising to 40 degrees once all satellites in
# * this shell are launched and operational.
# * 
# * Planned to use intersatellite links, which are not yet simulated
# * here - but the test duo and first batches of launched satellites
# * do not include those links.
#
# http://en.wikipedia.org/wiki/Starlink_(satellite_constellation)
#
# $Id$

# November 2018
# This is Orbital Plane 2 detailed in SAT-MOD-2018-1108-00083, FCC
# filing of 8 November 2018, and described in Appendix A.
# set SATS_PER_PLANE 66
# set NUM_PLANES 24

# August 2019
# This is detailed in the Technical Attachment of
# FCC document SAT-MOD-20190830-00087, filed 30 August 2019.
# SpaceX wants to rearrange its Starlink satellites for faster
# broadband ramp-up, Alan Boyle, GeekWire, 16 September 2019.
# https://www.geekwire.com/2019/spacex-wants-rearrange-starlink-satellites-faster-broadband-ramp/

set SATS_PER_PLANE 22
set NUM_PLANES 72

# to see planned initial deployment of third of shell,
# which provides continuous coverage at highest
# degrees of latitude, set this to 1
set DEPLOY_THIRD 0

# setup orbital elements
set a [expr 540.0 + $RADIUS_OF_EARTH]
set e 0.0
set inc 53.2
set omega 0.0
set T_per [expr 2 * $PI * pow($a,1.5) / sqrt($MU)]

# Minimum elevation angle is said to start at 25 degrees, rising to 40 degrees.
# Let's assume that 25 degrees is during partial constellation deployment,
# and 40 degrees is for this fully deployed constellation. Given uneven
# coverage for 40 degrees for the specified geometry, this seems doubtful.
# There are indications that, without intersatellite links, the lower
# minimum elevation of 25 degrees will be in use permanently.
set coverage_angle 25.0

upvar #0 NUM_COLORS NUM_COLORS

if {$NUM_COLORS < 19} {
    # more than 19 satellites can be visible in fisheye over mid-latitudes
    puts stderr "\nSaVi: Coverage view of Starlink constellation benefits from largest number of colors (19+)."
}

# Plane offset is really a function of harmonic factor in Ballard constellations.
# (Rosette Constellations of Earth Satellites, A. H. Ballard, TRW,
# IEEE Transactions on Aerospace and Electronic Systems, Vol 16 No 5, Sep. 1980)
# 360 / 66 / 24 = 360 / 1,584 = 0.2273 degrees approx.

# The Access .ndb database file in the November FCC filing suggests a
# ~3 degree offset between planes - 13th harmonic is 2.9545, close enough.
# set interplane_phasing [expr 360.0 / $NUM_PLANES / $SATS_PER_PLANE * 13]

# With the August 2019 filing, phasing is not specified.
# In the spirit of Starlink's constellation design to date, we'll just
# guess. Let's pick a prime and hope for the best. What could
# possibly go wrong? Offset makes so much more sense when negativa,
# since it's time TO periapsis.
set interplane_phasing [expr 360.0 / $NUM_PLANES / $SATS_PER_PLANE * -13]

satellites GV_BEGIN
for {set j 0} {$j < $NUM_PLANES} {incr j} {
	set Omega [expr $j * 360.0 / $NUM_PLANES]
	if {$DEPLOY_THIRD && [expr $j % 3 ] != 0} {
	          continue; 
	}
	for {set i 0} {$i < $SATS_PER_PLANE} {incr i} {
		set plane_offset [expr ($T_per / 360) * ($j * $interplane_phasing) ]

		set T [expr ($T_per * $i / $SATS_PER_PLANE ) + $plane_offset ]
		set n [satellites LOAD $a $e $inc $Omega $omega $T "Starlink-phase1-540km ($j, $i)"]
		if {$i > 0} {satellites ORBIT_SET $n 0}
	}
}
satellites GV_END

/// \file collision.c
/// \brief Predict and perform collisions between particles.
//
//	Copyright (c) 2014, Christian B. Mendl
//	All rights reserved.
//	http://christian.mendl.net
//
//	This program is free software; you can redistribute it and/or
//	modify it under the terms of the Simplified BSD License
//	http://www.opensource.org/licenses/bsd-license.php
//
//	References:
//	- Christian B. Mendl, Herbert Spohn
//	  Current fluctuations for anharmonic chains in thermal equilibrium
//	  arXiv:1412.4609
//
//	- Christian B. Mendl, Herbert Spohn
//	  Equilibrium time-correlation functions for one-dimensional hard-point systems
//	  Phys. Rev. E 90, 012147 (2014), arXiv:1403.0213
//
//	- Christian B. Mendl, Herbert Spohn
//	  Dynamic correlators of Fermi-Pasta-Ulam chains and nonlinear fluctuating hydrodynamics
//	  Phys. Rev. Lett. 111, 230601 (2013), arXiv:1305.1209
//
//	- Herbert Spohn
//	  Nonlinear fluctuating hydrodynamics for anharmonic chains
//	  J. Stat. Phys. 154, 1191-1227 (2014), arXiv:1305.6412
//_______________________________________________________________________________________________________________________
//

#include "collision.h"
#include "potential.h"
#include "util.h"
#include <math.h>
#include <stdlib.h>
#include <assert.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Predict next event (i.e., collision) between particles 'iprt' and 'iprt+1'
///
void PredictEvent(const field_t *f, const int iprt, event_t *ev)
{
	// initially inactivate event
	ev->active = false;

	// current lattice point
	const lattpoint_t *point = &f->pts[iprt];

	// time derivative of distance between neighboring particles, i.e., velocity
	// d r_i / dt = p_{i+1} - p_i
	double v = f->pts[(iprt + 1) & SITE_MASK].p - point->p;

	if (point->r > V_WIDTH_PLAT)	// particles are outside shoulder potential
	{
		if (v < 0)	// particles are approaching each other
		{
			// check whether kinetic energy is sufficient to enter shoulder plateau
			ev->type = square(v) - 4*V_POT_HEIGHT >= 0 ? COLL_PLAT_IN : COLL_PLAT_BOUNCE;
			ev->iprt = iprt;
			// calculate time point of collision
			ev->time = point->t + (V_WIDTH_PLAT - point->r) / v;
			ev->active = true;
		}
	}
	else	// point->r <= V_WIDTH_PLAT, particles are inside shoulder potential
	{
		assert(point->r >= V_WIDTH_CORE);

		if (v < 0)	// particles are approaching each other
		{
			ev->type = COLL_CORE;
			ev->iprt = iprt;
			// calculate time point of hard core collision
			ev->time = point->t + (V_WIDTH_CORE - point->r)/v;
			ev->active = true;
		}
		else if (v > 0)		// particles are separating
		{
			ev->type = COLL_PLAT_OUT;
			ev->iprt = iprt;
			// calculate time point of soft collision when leaving shoulder plateau
			ev->time = point->t + (V_WIDTH_PLAT - point->r)/v;
			ev->active = true;
		}
		else
		{
			// v == 0 should not occur
			assert(false);
		}
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Perform collision event 'ev' and update field variables
///
current_coll_t PerformCollision(const event_t *ev, field_t *f)
{
	assert(ev->active);

	// insert event into event log
	#ifdef _DEBUG
	EventLog_AddEvent(ev, &event_log);
	#endif

	// particle index
	int i = ev->iprt;

	// current lattice point
	lattpoint_t *point = &f->pts[i];

	// next and previous particles
	lattpoint_t *pt_next = &f->pts[(i + 1) & SITE_MASK];
	lattpoint_t *pt_prev = &f->pts[(i - 1) & SITE_MASK];

	// event must be in the future
	assert(ev->time >= point->t);

	// elongation current is negative time-integrated momentum up to 'ev->time'
	current_coll_t J;
	J.r      = (  point->t - ev->time) * point->p;
	J.r_next = (pt_next->t - ev->time) * pt_next->p;
	J.r_prev = (pt_prev->t - ev->time) * pt_prev->p;

	// advance positions of neighbors up to 'ev->time' since relevant momenta will change
	pt_next->r += (ev->time - pt_next->t) * (f->pts[(i+2) & SITE_MASK].p - pt_next->p);
	pt_prev->r += (ev->time - pt_prev->t) * (                   point->p - pt_prev->p);
	pt_next->t = ev->time;
	pt_prev->t = ev->time;

	// time derivative of distance between neighboring particles, i.e., velocity
	// d r_i / dt = p_{i+1} - p_i
	double v = pt_next->p - point->p;

	// store momenta for current calculation below
	double p0      = point->p;
	double p0_next = pt_next->p;

	if (ev->type == COLL_PLAT_IN)		// soft collision of approaching particles entering the shoulder plateau
	{
		// consistency checks:
		// particles must be approaching each other
		assert(v < 0);
		// current distance of particles must be larger than shoulder width
		assert(point->r > V_WIDTH_PLAT);
		// time point of collision
		assert(fabs((V_WIDTH_PLAT - point->r) - (ev->time - point->t)*v) < 1e-12);
		// kinetic energy must be sufficient to enter shoulder plateau
		assert(square(v) - 4*V_POT_HEIGHT >= 0);

		// at 'ev->time', particle distance is 'V_WIDTH_PLAT'
		point->r = V_WIDTH_PLAT;
		point->t = ev->time;

		// update momenta after collision
		double p_mean  = 0.5 * (pt_next->p + point->p);
		double p_delta = 0.5 * sqrt(square(v) - 4*V_POT_HEIGHT);
		point->p   = p_mean + p_delta;		// particles are approaching even after collision
		pt_next->p = p_mean - p_delta;
	}
	else if (ev->type == COLL_PLAT_OUT)		// soft collision of separating particles leaving the shoulder plateau
	{
		// consistency checks:
		// particles must be separating
		assert(v > 0);
		// current distance of particles must be within shoulder plateau
		assert(V_WIDTH_CORE <= point->r && point->r <= V_WIDTH_PLAT);
		// time point of collision
		assert(fabs((V_WIDTH_PLAT - point->r) - (ev->time - point->t)*v) < 1e-12);

		// at 'ev->time', particle distance is 'V_WIDTH_PLAT'
		point->r = V_WIDTH_PLAT;
		point->t = ev->time;

		// update momenta after collision
		double p_mean  = 0.5 * (pt_next->p + point->p);
		double p_delta = 0.5 * sqrt(square(v) + 4*V_POT_HEIGHT);
		point->p   = p_mean - p_delta;		// particles are still separating after collision
		pt_next->p = p_mean + p_delta;
	}
	else if (ev->type == COLL_PLAT_BOUNCE)	// bounce-back collision at shoulder plateau for approaching particles with small kinetic energy
	{
		// consistency checks:
		// particles must be approaching each other
		assert(v < 0);
		// current distance of particles must be larger than shoulder width
		assert(point->r > V_WIDTH_PLAT);
		// time point of collision
		assert(fabs((V_WIDTH_PLAT - point->r) - (ev->time - point->t)*v) < 1e-12);
		// kinetic energy must be insufficient to enter shoulder plateau
		assert(square(v) - 4*V_POT_HEIGHT < 0);

		// at 'ev->time', particle distance is 'V_WIDTH_PLAT'
		point->r = V_WIDTH_PLAT;
		point->t = ev->time;

		// simply swap momenta since masses are equal
		double curp = point->p;
		point->p    = pt_next->p;
		pt_next->p  = curp;
	}
	else if (ev->type == COLL_CORE)		// hard collision at core
	{
		// consistency checks:
		// particles must be approaching each other
		assert(v < 0);
		// current distance of particles must be within shoulder plateau
		assert(V_WIDTH_CORE <= point->r && point->r <= V_WIDTH_PLAT);
		// time point of collision
		assert(fabs((V_WIDTH_CORE - point->r) - (ev->time - point->t)*v) < 1e-12);

		// at 'ev->time', particle distance is 'V_WIDTH_CORE'
		point->r = V_WIDTH_CORE;
		point->t = ev->time;

		// simply swap momenta since masses are equal
		double curp = point->p;
		point->p    = pt_next->p;
		pt_next->p  = curp;
	}
	else
	{
		// unknown collision type
		assert(false);
	}

	// advance particle position by a small time step such that event at ev->time is not predicted again
	const double delta_t = 1e-3;
	point->r += delta_t * (pt_next->p - point->p);
	J.r -= delta_t * point->p;
	point->t += delta_t;

	// must invalidate events associated with neighbors since momenta have changed, and predict new events

	// calculate momentum and energy current
	J.p_next = 0.5*((pt_next->p - point->p) - (p0_next - p0));
	J.e_next = 0.5*(p0_next + p0)*J.p_next;
	if (ev->type == COLL_PLAT_IN)
	{
		J.e_next -= 0.5*V_POT_HEIGHT;
	}
	else if (ev->type == COLL_PLAT_OUT)
	{
		J.e_next += 0.5*V_POT_HEIGHT;
	}
	return J;
}



#ifdef _DEBUG

//_______________________________________________________________________________________________________________________
///
/// \brief Insert event into event log
///
void EventLog_AddEvent(const event_t *ev, event_log_t *evlog)
{
	if (evlog->last != NULL)
	{
		assert(evlog->num > 0 && evlog->first != NULL);
		evlog->last->next = malloc(sizeof(event_node_t));
		evlog->last = evlog->last->next;
	}
	else
	{
		// event log is still empty
		assert(evlog->num == 0 && evlog->first == NULL);
		evlog->first = evlog->last = malloc(sizeof(event_node_t));
	}

	// copy event
	evlog->last->ev = *ev;
	evlog->last->ev.time += evlog->time_offset;	// add time offset
	evlog->last->next = NULL;

	evlog->num++;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Delete event log (free memory)
///
void EventLog_Delete(event_log_t *evlog)
{
	while (evlog->first != NULL)
	{
		event_node_t *next = evlog->first->next;
		free(evlog->first);
		evlog->first = next;
		evlog->num--;
	}
	assert(evlog->num == 0);

	evlog->last = NULL;
	evlog->time_offset = 0;
}

//_______________________________________________________________________________________________________________________
///
/// \brief Global event log for testing
///
event_log_t event_log;


#endif	// _DEBUG

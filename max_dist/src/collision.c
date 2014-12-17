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

	assert(0 <= point->r && point->r <= V_WIDTH_WELL);

	// relative velocity
	double u = f->pts[(iprt + 1) & SITE_MASK].v - point->v;

	if (u < 0)	// particles are approaching each other
	{
		ev->type = COLL_CORE;
		ev->iprt = iprt;
		// calculate time point of hard-core collision
		ev->time = point->t - point->r/u;
		ev->active = true;
	}
	else if (u > 0)		// particles are separating
	{
		ev->type = COLL_STRING;
		ev->iprt = iprt;
		// calculate time point of of "collision" at distance 'V_WIDTH_WELL'
		ev->time = point->t + (V_WIDTH_WELL - point->r)/u;
		ev->active = true;
	}
	else
	{
		// u == 0 should not occur
		assert(false);
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

	assert(0 <= point->r && point->r <= V_WIDTH_WELL);

	// next and previous particles
	lattpoint_t *pt_next = &f->pts[(i + 1) & SITE_MASK];
	lattpoint_t *pt_prev = &f->pts[(i - 1) & SITE_MASK];

	// event must be in the future
	assert(ev->time >= point->t);

	// elongation current is negative time-integrated velocity up to 'ev->time'
	current_coll_t J;
	J.r      = (  point->t - ev->time) * point->v;
	J.r_next = (pt_next->t - ev->time) * pt_next->v;
	J.r_prev = (pt_prev->t - ev->time) * pt_prev->v;

	// advance positions of neighbors up to 'ev->time' since relevant momenta will change
	pt_next->r += (ev->time - pt_next->t) * (f->pts[(i+2) & SITE_MASK].v - pt_next->v);
	pt_prev->r += (ev->time - pt_prev->t) * (                   point->v - pt_prev->v);
	pt_next->t = ev->time;
	pt_prev->t = ev->time;

	// relative velocity
	double u = pt_next->v - point->v;

	if (ev->type == COLL_CORE)		// elastic collision of approaching particles at core
	{
		// consistency checks:
		// particles must be approaching each other
		assert(u < 0);
		// time point of collision
		assert(fabs(point->r + (ev->time - point->t)*u) < 1e-12);

		// at 'ev->time', particle distance is zero
		point->r = 0;
		point->t = ev->time;
	}
	else if (ev->type == COLL_STRING)		// elastic "collision" at distance 'V_WIDTH_WELL' effected by virtual string
	{
		// consistency checks:
		// particles must be separating
		assert(u > 0);
		// time point of collision
		assert(fabs((V_WIDTH_WELL - point->r) - (ev->time - point->t)*u) < 1e-12);

		// at 'ev->time', particle distance is 'V_WIDTH_WELL'
		point->r = V_WIDTH_WELL;
		point->t = ev->time;
	}
	else
	{
		// unknown collision type
		assert(false);
	}

	// total velocity of system consisting of current and next point
	double v_tot;

	// elastic collision
	if (i % 2 == 0)		// current point has mass m_0
	{
		v_tot = 0.5*(f->relm[0]*point->v + f->relm[1]*pt_next->v);

		point->v   += f->relm[1]*u;
		pt_next->v -= f->relm[0]*u;
	}
	else				// current point has mass m_1
	{
		v_tot = 0.5*(f->relm[1]*point->v + f->relm[0]*pt_next->v);

		point->v   += f->relm[0]*u;
		pt_next->v -= f->relm[1]*u;
	}

	// advance particle position by a small time step such that event at ev->time is not predicted again
	const double delta_t = 1e-3;
	point->r += delta_t * (pt_next->v - point->v);
	J.r -= delta_t * point->v;
	point->t += delta_t;

	// must invalidate events associated with neighbors since momenta have changed, and predict new events

	// calculate momentum and energy current
	J.p_next = -f->rm_prod*u;
	J.e_next = v_tot*J.p_next;
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

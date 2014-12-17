/// \file simulation.c
/// \brief "Sweep" the event table, i.e., process all active events within current time period.
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

#include "simulation.h"
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Collect all collision events for field 'f' (used to initialize event table 'tab')
///
void CollectEvents(field_t *f, event_table_t *tab)
{
	unsigned int iprt;
	for (iprt = 0; iprt < NUM_SITES; iprt++)
	{
		event_t ev;
		PredictEvent(f, iprt, &ev);
		if (ev.active)
		{
			assert(ev.time >= 0);

			// insert event into table and store event index
			f->pts[iprt].ievent = EventTable_InsertEvent(&ev, tab);
		}
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Update event associated with particle 'iprt' (remove associated event and predict new event)
///
static void UpdateFieldEvent(const int iprt, field_t *f, event_table_t *tab)
{
	// remove event associated with particle, if any
	if (f->pts[iprt].ievent >= 0) {
		EventTable_RemoveEvent(f->pts[iprt].ievent, tab);
		// will be set below
		// f->pts[iprt].ievent = -1;
	}

	// predict new event, if any
	event_t ev;
	PredictEvent(f, iprt, &ev);
	f->pts[iprt].ievent = EventTable_InsertEvent(&ev, tab);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Process all active events within current time period and return total current
///
current_t SweepTable(field_t *f, event_table_t *tab)
{
	// time period
	const double period = NUM_SLOTS * tab->dt;

	// record cumulative momentum and energy current
	current_t J = { 0 };

	unsigned int islot;
	for (islot = 0; islot < NUM_SLOTS; islot++)
	{
		// reference to time slot
		event_slot_t *slot = &tab->slots[islot];

		while (slot->num > 0)	// as long as there are active events in current time slot
		{
			// process and update all events within current period

			// find earliest active event within time slot
			double tmin = period;	// ignore events beyond current time period
			int imin = -1;
			int i;
			for (i = 0; i < EVENTS_PER_SLOT; i++)
			{
				if (slot->events[i].active && slot->events[i].time < tmin)
				{
					tmin = slot->events[i].time;
					imin = i;
				}
			}
			// active event found?
			if (imin < 0) {
				break;	// next slot
			}

			// reference to event
			const event_t *ev = &slot->events[imin];

			// consistency check: particle of event must link back to event
			assert(f->pts[ev->iprt].ievent == (int)islot*EVENTS_PER_SLOT + imin);
			// consistency check: appropriate time slot
			assert(islot*tab->dt <= ev->time && ev->time < (islot + 1)*tab->dt);

			current_coll_t J_cur = PerformCollision(ev, f);
			// accumulate current
			J.r += J_cur.r_prev + J_cur.r + J_cur.r_next;
			J.p += J_cur.p_next;
			J.e += J_cur.e_next;

			// copy particle index since event pointer might be overwritten by 'UpdateFieldEvent' below
			int iprt = ev->iprt;

			// update events associated with current particle and neighbors since momenta might have changed;
			// note that this process might re-insert events into current time slot
			UpdateFieldEvent( iprt,                  f, tab);
			UpdateFieldEvent((iprt + 1) & SITE_MASK, f, tab);
			UpdateFieldEvent((iprt - 1) & SITE_MASK, f, tab);
		}
	}

	// normalize total current
	J.r /= (NUM_SITES * period);
	J.p /= (NUM_SITES * period);
	J.e /= (NUM_SITES * period);

	return J;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Process all active events within current time period and return current for each lattice site
///
void SweepTableJSite(field_t *f, event_table_t *tab, current_t J[NUM_SITES])
{
	// time period
	const double period = NUM_SLOTS * tab->dt;

	// initially set to zero
	memset(J, 0, NUM_SITES * sizeof(current_t));

	unsigned int islot;
	for (islot = 0; islot < NUM_SLOTS; islot++)
	{
		// reference to time slot
		event_slot_t *slot = &tab->slots[islot];

		while (slot->num > 0)	// as long as there are active events in current time slot
		{
			// process and update all events within current period

			// find earliest active event within time slot
			double tmin = period;	// ignore events beyond current time period
			int imin = -1;
			int i;
			for (i = 0; i < EVENTS_PER_SLOT; i++)
			{
				if (slot->events[i].active && slot->events[i].time < tmin)
				{
					tmin = slot->events[i].time;
					imin = i;
				}
			}
			// active event found?
			if (imin < 0) {
				break;	// next slot
			}

			// reference to event
			const event_t *ev = &slot->events[imin];

			// consistency check: particle of event must link back to event
			assert(f->pts[ev->iprt].ievent == (int)islot*EVENTS_PER_SLOT + imin);
			// consistency check: appropriate time slot
			assert(islot*tab->dt <= ev->time && ev->time < (islot + 1)*tab->dt);

			// copy particle index since event pointer might be overwritten by 'UpdateFieldEvent' below
			int iprt = ev->iprt;
			int iprt_next = (iprt + 1) & SITE_MASK;
			int iprt_prev = (iprt - 1) & SITE_MASK;

			current_coll_t J_cur = PerformCollision(ev, f);
			// accumulate current
			J[iprt_prev].r += J_cur.r_prev;
			J[iprt     ].r += J_cur.r;
			J[iprt_next].r += J_cur.r_next;
			J[iprt_next].p += J_cur.p_next;
			J[iprt_next].e += J_cur.e_next;

			// update events associated with current particle and neighbors since momenta might have changed;
			// note that this process might re-insert events into current time slot
			UpdateFieldEvent( iprt,                  f, tab);
			UpdateFieldEvent((iprt + 1) & SITE_MASK, f, tab);
			UpdateFieldEvent((iprt - 1) & SITE_MASK, f, tab);
		}
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Subtract 't_shift' from all timing events (to avoid overflow)
///
void GlobalTimeshift(const double t_shift, field_t *f, event_table_t *tab)
{
	// event table
	unsigned int islot;
	for (islot = 0; islot < NUM_SLOTS; islot++)
	{
		// reference to time slot
		event_slot_t *slot = &tab->slots[islot];

		int i;
		for (i = 0; i < EVENTS_PER_SLOT; i++) {
			slot->events[i].time -= t_shift;
		}
	}

	// lattice field
	unsigned int iprt;
	for (iprt = 0; iprt < NUM_SITES; iprt++)
	{
		f->pts[iprt].t -= t_shift;
	}
}

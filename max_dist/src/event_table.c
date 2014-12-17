/// \file event_table.c
/// \brief Event table for collecting collision events in time slots.
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

#include "event_table.h"
#include <fftw3.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <stdio.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Allocate memory for event table
///
void EventTable_Allocate(const double dt, event_table_t *tab)
{
	tab->dt = dt;
	// aligned memory allocation
	tab->slots = fftw_malloc(NUM_SLOTS * sizeof(event_slot_t));
	memset(tab->slots, 0, NUM_SLOTS * sizeof(event_slot_t));	// all events are initially inactive
	tab->numevents = 0;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Delete event table (free memory)
///
void EventTable_Delete(event_table_t *tab)
{
	tab->numevents = 0;
	fftw_free(tab->slots);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Reset event table
///
void EventTable_Reset(event_table_t *tab)
{
	memset(tab->slots, 0, NUM_SLOTS * sizeof(event_slot_t));	// make all events inactive
	tab->numevents = 0;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Insert an event into table
///
/// \param ev event to be inserted
/// \param tab event table
/// \returns index of the inserted event, or -1 if event was not inserted
///
int EventTable_InsertEvent(const event_t *ev, event_table_t *tab)
{
	int i;

	// inactive events are ignored
	if (!ev->active) {
		return -1;
	}

	assert(ev->time >= 0);

	// find appropriate time slot;
	// 'unsigned int' to avoid invalid negative index after integer overflow
	unsigned int islot = (unsigned int)(ev->time / tab->dt) & SLOT_MASK;	// cast to integer to get floor value; cyclic index

	// reference to time slot
	event_slot_t *slot = &tab->slots[islot];

	// find an inactive event which can be overwritten
	for (i = 0; i < EVENTS_PER_SLOT; i++)
	{
		if (!slot->events[i].active)
		{
			// overwrite event
			slot->events[i] = *ev;
			// increase event counter
			slot->num++;
			tab->numevents++;
			// return event index
			return islot*EVENTS_PER_SLOT + i;
		}
	}

	// event slot overflow, not implemented yet
	fprintf(stderr, "event slot overflow, maximum number of events per slot: %i\n", EVENTS_PER_SLOT);
	assert(slot->num == EVENTS_PER_SLOT);
	exit(-1);

	//return -1;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Remove event from table
///
void EventTable_RemoveEvent(const int ievent, event_table_t *tab)
{
	assert(ievent >= 0);

	unsigned int islot = ievent / EVENTS_PER_SLOT;
	assert(islot < NUM_SLOTS);

	// reference to time slot
	event_slot_t *slot = &tab->slots[islot];

	// event should be active
	assert(slot->events[ievent & EVENT_MASK].active);

	// inactivate event
	slot->events[ievent & EVENT_MASK].active = false;

	// decrease event count
	slot->num--;
	tab->numevents--;
}

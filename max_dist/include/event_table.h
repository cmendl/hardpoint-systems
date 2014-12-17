//	Event table for collecting collision events in time slots.
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

#ifndef EVENT_TABLE_H
#define EVENT_TABLE_H

#include "field.h"
#include "collision.h"


//_______________________________________________________________________________________________________________________
//

#define NUM_SLOTS			8192					//!< number of event slots
#define SLOT_MASK			(NUM_SLOTS - 1)			//!< bit mask for fast modulo NUM_SLOTS operation; NUM_SLOTS must be a power of 2

#define EVENTS_PER_SLOT		16						//!< number of events per slot
#define EVENT_MASK			(EVENTS_PER_SLOT - 1)	//!< bit mask for fast modulo EVENTS_PER_SLOT operation; EVENTS_PER_SLOT must be a power of 2


//_______________________________________________________________________________________________________________________
///
/// \brief Time slot for storing events
///
typedef struct
{
	event_t events[EVENTS_PER_SLOT];	//!< events in the slot
	unsigned int num;					//!< number of active events
}
event_slot_t;


//_______________________________________________________________________________________________________________________
///
/// \brief Event table for collecting collision events in time slots
///
/// After a complete sweep of the table, the full time interval 'nslots * dt'
/// is subtracted from all timing events to avoid overflow.
///
typedef struct
{
	double dt;				//!< length (time interval) of a time slot
	event_slot_t *slots;	//!< event slots, time points are interpreted modulo 'nslots * dt'
	unsigned int numevents;	//!< total number of events in table
}
event_table_t;


void EventTable_Allocate(const double dt, event_table_t *tab);

void EventTable_Delete(event_table_t *tab);

void EventTable_Reset(event_table_t *tab);


//_______________________________________________________________________________________________________________________
//

int EventTable_InsertEvent(const event_t *ev, event_table_t *tab);

void EventTable_RemoveEvent(const int ievent, event_table_t *tab);



#endif

//	Predict and perform collisions between particles.
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

#ifndef COLLISION_H
#define COLLISION_H

#include "field.h"
#include "current.h"
#include <stdbool.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Collision type
///
typedef enum
{
	COLL_STRING,		//!< elastic "collision" at distance 'V_WIDTH_WELL' effected by virtual string
	COLL_CORE,			//!< elastic collision of approaching particles at core
}
colltype_t;


//_______________________________________________________________________________________________________________________
///
/// \brief Collision event
///
typedef struct
{
	double time;		//!< time point of the event
	int iprt;			//!< particle index
	colltype_t type;	//!< collision type
	bool active;		//!< whether the event is active
}
event_t;


void PredictEvent(const field_t *f, const int iprt, event_t *ev);


//_______________________________________________________________________________________________________________________
///
/// \brief current J of various lattice fields during collision
///
typedef struct
{
	double r_prev;		//!< elongation current of previous site
	double r;			//!< elongation current of actual site
	double r_next;		//!< elongation current of next site
	double p_next;		//!< momentum current of next site
	double e_next;		//!< energy current of next site
}
current_coll_t;


//_______________________________________________________________________________________________________________________
//


current_coll_t PerformCollision(const event_t *ev, field_t *f);



#ifdef _DEBUG

//_______________________________________________________________________________________________________________________
///
/// \brief Node of linked list for event log
///
typedef struct event_node_s
{
	event_t ev;
	struct event_node_s *next;
}
event_node_t;

//_______________________________________________________________________________________________________________________
///
/// \brief Event log for testing (linked list)
///
typedef struct
{
	event_node_t *first;	//!< pointer to first event, or NULL if there aren't any events
	event_node_t *last;		//!< pointer to last event
	unsigned int num;		//!< number of events
	double time_offset;		//!< time offset to compensate the global timeshift after sweeps of the event table
}
event_log_t;

void EventLog_AddEvent(const event_t *ev, event_log_t *evlog);

void EventLog_Delete(event_log_t *evlog);


// global event log for testing
extern event_log_t event_log;


#endif	// _DEBUG


#endif

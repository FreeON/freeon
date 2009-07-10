/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* Programmer:  John Mainzer
 *              10/27/05
 *
 *		This file contains common code for tests of the cache
 *		implemented in H5C.c
 */
#include "h5test.h"
#include "H5Iprivate.h"
#include "H5ACprivate.h"
#include "cache_common.h"


/* global variable declarations: */

hbool_t write_permitted = TRUE;
hbool_t pass = TRUE; /* set to false on error */
hbool_t skip_long_tests = TRUE;
hbool_t run_full_test = TRUE;
const char *failure_mssg = NULL;

test_entry_t pico_entries[NUM_PICO_ENTRIES];
test_entry_t nano_entries[NUM_NANO_ENTRIES];
test_entry_t micro_entries[NUM_MICRO_ENTRIES];
test_entry_t tiny_entries[NUM_TINY_ENTRIES];
test_entry_t small_entries[NUM_SMALL_ENTRIES];
test_entry_t medium_entries[NUM_MEDIUM_ENTRIES];
test_entry_t large_entries[NUM_LARGE_ENTRIES];
test_entry_t huge_entries[NUM_HUGE_ENTRIES];
test_entry_t monster_entries[NUM_MONSTER_ENTRIES];
test_entry_t variable_entries[NUM_VARIABLE_ENTRIES];

test_entry_t * entries[NUMBER_OF_ENTRY_TYPES] =
{
    pico_entries,
    nano_entries,
    micro_entries,
    tiny_entries,
    small_entries,
    medium_entries,
    large_entries,
    huge_entries,
    monster_entries,
    variable_entries
};

const int32_t max_indices[NUMBER_OF_ENTRY_TYPES] =
{
    NUM_PICO_ENTRIES - 1,
    NUM_NANO_ENTRIES - 1,
    NUM_MICRO_ENTRIES - 1,
    NUM_TINY_ENTRIES - 1,
    NUM_SMALL_ENTRIES - 1,
    NUM_MEDIUM_ENTRIES - 1,
    NUM_LARGE_ENTRIES - 1,
    NUM_HUGE_ENTRIES - 1,
    NUM_MONSTER_ENTRIES - 1,
    NUM_VARIABLE_ENTRIES - 1
};

const size_t entry_sizes[NUMBER_OF_ENTRY_TYPES] =
{
    PICO_ENTRY_SIZE,
    NANO_ENTRY_SIZE,
    MICRO_ENTRY_SIZE,
    TINY_ENTRY_SIZE,
    SMALL_ENTRY_SIZE,
    MEDIUM_ENTRY_SIZE,
    LARGE_ENTRY_SIZE,
    HUGE_ENTRY_SIZE,
    MONSTER_ENTRY_SIZE,
    VARIABLE_ENTRY_SIZE
};

const haddr_t base_addrs[NUMBER_OF_ENTRY_TYPES] =
{
    PICO_BASE_ADDR,
    NANO_BASE_ADDR,
    MICRO_BASE_ADDR,
    TINY_BASE_ADDR,
    SMALL_BASE_ADDR,
    MEDIUM_BASE_ADDR,
    LARGE_BASE_ADDR,
    HUGE_BASE_ADDR,
    MONSTER_BASE_ADDR,
    VARIABLE_BASE_ADDR
};

const haddr_t alt_base_addrs[NUMBER_OF_ENTRY_TYPES] =
{
    PICO_ALT_BASE_ADDR,
    NANO_ALT_BASE_ADDR,
    MICRO_ALT_BASE_ADDR,
    TINY_ALT_BASE_ADDR,
    SMALL_ALT_BASE_ADDR,
    MEDIUM_ALT_BASE_ADDR,
    LARGE_ALT_BASE_ADDR,
    HUGE_ALT_BASE_ADDR,
    MONSTER_ALT_BASE_ADDR,
    VARIABLE_ALT_BASE_ADDR
};

const char * entry_type_names[NUMBER_OF_ENTRY_TYPES] =
{
    "pico entries -- 1 B",
    "nano entries -- 4 B",
    "micro entries -- 16 B",
    "tiny entries -- 64 B",
    "small entries -- 256 B",
    "medium entries -- 1 KB",
    "large entries -- 4 KB",
    "huge entries -- 16 KB",
    "monster entries -- 64 KB",
    "variable entries -- 1B - 10KB"
};


/* callback table declaration */

const H5C_class_t types[NUMBER_OF_ENTRY_TYPES] =
{
  {
    PICO_ENTRY_TYPE,
    (H5C_load_func_t)pico_load,
    (H5C_flush_func_t)pico_flush,
    (H5C_dest_func_t)pico_dest,
    (H5C_clear_func_t)pico_clear,
    (H5C_size_func_t)pico_size
  },
  {
    NANO_ENTRY_TYPE,
    (H5C_load_func_t)nano_load,
    (H5C_flush_func_t)nano_flush,
    (H5C_dest_func_t)nano_dest,
    (H5C_clear_func_t)nano_clear,
    (H5C_size_func_t)nano_size
  },
  {
    MICRO_ENTRY_TYPE,
    (H5C_load_func_t)micro_load,
    (H5C_flush_func_t)micro_flush,
    (H5C_dest_func_t)micro_dest,
    (H5C_clear_func_t)micro_clear,
    (H5C_size_func_t)micro_size
  },
  {
    TINY_ENTRY_TYPE,
    (H5C_load_func_t)tiny_load,
    (H5C_flush_func_t)tiny_flush,
    (H5C_dest_func_t)tiny_dest,
    (H5C_clear_func_t)tiny_clear,
    (H5C_size_func_t)tiny_size
  },
  {
    SMALL_ENTRY_TYPE,
    (H5C_load_func_t)small_load,
    (H5C_flush_func_t)small_flush,
    (H5C_dest_func_t)small_dest,
    (H5C_clear_func_t)small_clear,
    (H5C_size_func_t)small_size
  },
  {
    MEDIUM_ENTRY_TYPE,
    (H5C_load_func_t)medium_load,
    (H5C_flush_func_t)medium_flush,
    (H5C_dest_func_t)medium_dest,
    (H5C_clear_func_t)medium_clear,
    (H5C_size_func_t)medium_size
  },
  {
    LARGE_ENTRY_TYPE,
    (H5C_load_func_t)large_load,
    (H5C_flush_func_t)large_flush,
    (H5C_dest_func_t)large_dest,
    (H5C_clear_func_t)large_clear,
    (H5C_size_func_t)large_size
  },
  {
    HUGE_ENTRY_TYPE,
    (H5C_load_func_t)huge_load,
    (H5C_flush_func_t)huge_flush,
    (H5C_dest_func_t)huge_dest,
    (H5C_clear_func_t)huge_clear,
    (H5C_size_func_t)huge_size
  },
  {
    MONSTER_ENTRY_TYPE,
    (H5C_load_func_t)monster_load,
    (H5C_flush_func_t)monster_flush,
    (H5C_dest_func_t)monster_dest,
    (H5C_clear_func_t)monster_clear,
    (H5C_size_func_t)monster_size
  },
  {
    VARIABLE_ENTRY_TYPE,
    (H5C_load_func_t)variable_load,
    (H5C_flush_func_t)variable_flush,
    (H5C_dest_func_t)variable_dest,
    (H5C_clear_func_t)variable_clear,
    (H5C_size_func_t)variable_size
  }
};

static herr_t clear(H5F_t * f, void * thing, hbool_t dest);
static herr_t destroy(H5F_t * f, void * thing);
static herr_t flush(H5F_t *f, hid_t dxpl_id, hbool_t dest,
                    haddr_t addr, void *thing, unsigned UNUSED * flags_ptr);
static void * load(H5F_t *f, hid_t dxpl_id, haddr_t addr,
                   const void *udata1, void *udata2);
static herr_t size(H5F_t * f, void * thing, size_t * size_ptr);



/* address translation funtions: */

/*-------------------------------------------------------------------------
 * Function:	addr_to_type_and_index
 *
 * Purpose:	Given an address, compute the type and index of the
 *		associated entry.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              6/10/04
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
void
addr_to_type_and_index(haddr_t addr,
                       int32_t * type_ptr,
                       int32_t * index_ptr)
{
    int i;
    int32_t type;
    int32_t idx;

    HDassert( type_ptr );
    HDassert( index_ptr );

    /* we only have a small number of entry types, so just do a
     * linear search.  If NUMBER_OF_ENTRY_TYPES grows, we may want
     * to do a binary search instead.
     */
    i = 1;
    if ( addr >= PICO_ALT_BASE_ADDR ) {

        while ( ( i < NUMBER_OF_ENTRY_TYPES ) &&
                ( addr >= alt_base_addrs[i] ) )
        {
            i++;
        }

    } else {

        while ( ( i < NUMBER_OF_ENTRY_TYPES ) &&
                ( addr >= base_addrs[i] ) )
        {
            i++;
        }
    }

    type = i - 1;

    HDassert( ( type >= 0 ) && ( type < NUMBER_OF_ENTRY_TYPES ) );

    if ( addr >= PICO_ALT_BASE_ADDR ) {

        idx = (int32_t)((addr - alt_base_addrs[type]) / entry_sizes[type]);
        HDassert( ( idx >= 0 ) && ( idx <= max_indices[type] ) );
        HDassert( !((entries[type])[idx].at_main_addr) );
        HDassert( addr == (entries[type])[idx].alt_addr );

    } else {

        idx = (int32_t)((addr - base_addrs[type]) / entry_sizes[type]);
        HDassert( ( idx >= 0 ) && ( idx <= max_indices[type] ) );
        HDassert( (entries[type])[idx].at_main_addr );
        HDassert( addr == (entries[type])[idx].main_addr );
    }

    HDassert( addr == (entries[type])[idx].addr );

    *type_ptr = type;
    *index_ptr = idx;

    return;

} /* addr_to_type_and_index() */


#if 0 /* This function has never been used, but we may want it
       * some time.  Lets keep it for now.
       */
/*-------------------------------------------------------------------------
 * Function:	type_and_index_to_addr
 *
 * Purpose:	Given a type and index of an entry, compute the associated
 *		addr and return that value.
 *
 * Return:	computed addr
 *
 * Programmer:	John Mainzer
 *              6/10/04
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
haddr_t
type_and_index_to_addr(int32_t type,
                       int32_t idx)
{
    haddr_t addr;

    HDassert( ( type >= 0 ) && ( type < NUMBER_OF_ENTRY_TYPES ) );
    HDassert( ( idx >= 0 ) && ( idx <= max_indices[type] ) );

    addr = base_addrs[type] + (((haddr_t)idx) * entry_sizes[type]);

    HDassert( addr == (entries[type])[idx].addr );

    if ( (entries[type])[idx].at_main_addr ) {

        HDassert( addr == (entries[type])[idx].main_addr );

    } else {

        HDassert( addr == (entries[type])[idx].alt_addr );
    }

    return(addr);

} /* type_and_index_to_addr() */

#endif


/* Call back functions: */

/*-------------------------------------------------------------------------
 *
 * Function:    check_if_write_permitted
 *
 * Purpose:     Determine if a write is permitted under the current
 *              circumstances, and set *write_permitted_ptr accordingly.
 *              As a general rule it is, but when we are running in parallel
 *              mode with collective I/O, we must ensure that a read cannot
 *              cause a write.
 *
 *              In the event of failure, the value of *write_permitted_ptr
 *              is undefined.
 *
 * Return:      Non-negative on success/Negative on failure.
 *
 * Programmer:  John Mainzer, 5/15/04
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */

herr_t
check_write_permitted(const H5F_t UNUSED * f,
                      hid_t UNUSED dxpl_id,
                      hbool_t * write_permitted_ptr)
{

    HDassert( write_permitted_ptr );
    *write_permitted_ptr = write_permitted;

    return(SUCCEED);

} /* check_write_permitted() */


/*-------------------------------------------------------------------------
 * Function:	clear & friends
 *
 * Purpose:	clear the entry.  The helper functions verify that the
 *		correct version of clear is being called, and then call
 *		clear proper.
 *
 * Return:	SUCCEED
 *
 * Programmer:	John Mainzer
 *              6/10/04
 *
 * Modifications:
 *
 * 		Added variable_clear.  -- JRM 8/30/06
 *
 *-------------------------------------------------------------------------
 */

herr_t
clear(H5F_t * f,
      void *  thing,
      hbool_t dest)
{
    test_entry_t * entry_ptr;
    test_entry_t * base_addr;

    HDassert( thing );

    entry_ptr = (test_entry_t *)thing;
    base_addr = entries[entry_ptr->type];

    HDassert( entry_ptr->index >= 0 );
    HDassert( entry_ptr->index <= max_indices[entry_ptr->type] );
    HDassert( entry_ptr == &(base_addr[entry_ptr->index]) );
    HDassert( entry_ptr == entry_ptr->self );
    HDassert( entry_ptr->header.addr == entry_ptr->addr );
    HDassert( entry_ptr->header.size == entry_ptr->size );
    HDassert( ( entry_ptr->type == VARIABLE_ENTRY_TYPE ) ||
	      ( entry_ptr->size == entry_sizes[entry_ptr->type] ) );

    entry_ptr->header.is_dirty = FALSE;
    entry_ptr->is_dirty = FALSE;

    entry_ptr->cleared = TRUE;

    if ( dest ) {

        destroy(f, thing);

    }

    return(SUCCEED);

} /* clear() */

herr_t
pico_clear(H5F_t * f, void *  thing, hbool_t dest)
{
    HDassert ( ((test_entry_t *)thing)->type == PICO_ENTRY_TYPE );
    return(clear(f, thing, dest));
}

herr_t
nano_clear(H5F_t * f, void *  thing, hbool_t dest)
{
    HDassert ( ((test_entry_t *)thing)->type == NANO_ENTRY_TYPE );
    return(clear(f, thing, dest));
}

herr_t
micro_clear(H5F_t * f, void *  thing, hbool_t dest)
{
    HDassert ( ((test_entry_t *)thing)->type == MICRO_ENTRY_TYPE );
    return(clear(f, thing, dest));
}

herr_t
tiny_clear(H5F_t * f, void *  thing, hbool_t dest)
{
    HDassert ( ((test_entry_t *)thing)->type == TINY_ENTRY_TYPE );
    return(clear(f, thing, dest));
}

herr_t
small_clear(H5F_t * f, void *  thing, hbool_t dest)
{
    HDassert ( ((test_entry_t *)thing)->type == SMALL_ENTRY_TYPE );
    return(clear(f, thing, dest));
}

herr_t
medium_clear(H5F_t * f, void *  thing, hbool_t dest)
{
    HDassert ( ((test_entry_t *)thing)->type == MEDIUM_ENTRY_TYPE );
    return(clear(f, thing, dest));
}

herr_t
large_clear(H5F_t * f, void *  thing, hbool_t dest)
{
    HDassert ( ((test_entry_t *)thing)->type == LARGE_ENTRY_TYPE );
    return(clear(f, thing, dest));
}

herr_t
huge_clear(H5F_t * f, void *  thing, hbool_t dest)
{
    HDassert ( ((test_entry_t *)thing)->type == HUGE_ENTRY_TYPE );
    return(clear(f, thing, dest));
}

herr_t
monster_clear(H5F_t * f, void *  thing, hbool_t dest)
{
    HDassert ( ((test_entry_t *)thing)->type == MONSTER_ENTRY_TYPE );
    return(clear(f, thing, dest));
}

herr_t
variable_clear(H5F_t * f, void *  thing, hbool_t dest)
{
    HDassert ( ((test_entry_t *)thing)->type == VARIABLE_ENTRY_TYPE );
    return(clear(f, thing, dest));
}



/*-------------------------------------------------------------------------
 * Function:	dest & friends
 *
 * Purpose:	Destroy the entry.  The helper functions verify that the
 *		correct version of dest is being called, and then call
 *		dest proper.
 *
 * Return:	SUCCEED
 *
 * Programmer:	John Mainzer
 *              6/10/04
 *
 * Modifications:
 *
 * 		JRM -- 4/4/06
 * 		Added code to decrement the pinning_ref_count s of entries
 * 		pinned by the target entry, and to unpin those entries
 * 		if the reference count drops to zero.
 *
 * 		JRM -- 8/30/06
 * 		Added variable_destroy().
 *
 *-------------------------------------------------------------------------
 */

herr_t
destroy(H5F_t UNUSED * f,
        void *         thing)
{
    int i;
    test_entry_t * entry_ptr;
    test_entry_t * base_addr;
    test_entry_t * pinned_entry_ptr;
    test_entry_t * pinned_base_addr;

    HDassert( thing );

    entry_ptr = (test_entry_t *)thing;
    base_addr = entries[entry_ptr->type];

    HDassert( entry_ptr->index >= 0 );
    HDassert( entry_ptr->index <= max_indices[entry_ptr->type] );
    HDassert( entry_ptr == &(base_addr[entry_ptr->index]) );
    HDassert( entry_ptr == entry_ptr->self );
    HDassert( entry_ptr->cache_ptr != NULL );
    HDassert( entry_ptr->cache_ptr->magic == H5C__H5C_T_MAGIC );
    HDassert( ( entry_ptr->header.destroy_in_progress ) ||
              ( entry_ptr->header.addr == entry_ptr->addr ) );
    HDassert( entry_ptr->header.size == entry_ptr->size );
    HDassert( ( entry_ptr->type == VARIABLE_ENTRY_TYPE ) ||
	      ( entry_ptr->size == entry_sizes[entry_ptr->type] ) );

    HDassert( !(entry_ptr->is_dirty) );
    HDassert( !(entry_ptr->header.is_dirty) );

    if ( entry_ptr->num_pins > 0 ) {

	for ( i = 0; i < entry_ptr->num_pins; i++ )
        {
	    pinned_base_addr = entries[entry_ptr->pin_type[i]];
	    pinned_entry_ptr = &(pinned_base_addr[entry_ptr->pin_idx[i]]);

	    HDassert( 0 <= pinned_entry_ptr->type );
            HDassert( pinned_entry_ptr->type < NUMBER_OF_ENTRY_TYPES );
	    HDassert( pinned_entry_ptr->type == entry_ptr->pin_type[i] );
	    HDassert( pinned_entry_ptr->index >= 0 );
	    HDassert( pinned_entry_ptr->index <=
		      max_indices[pinned_entry_ptr->type] );
	    HDassert( pinned_entry_ptr->index == entry_ptr->pin_idx[i] );
	    HDassert( pinned_entry_ptr == pinned_entry_ptr->self );
	    HDassert( pinned_entry_ptr->header.is_pinned );
	    HDassert( pinned_entry_ptr->is_pinned );
	    HDassert( pinned_entry_ptr->pinning_ref_count > 0 );

	    pinned_entry_ptr->pinning_ref_count--;

	    if ( pinned_entry_ptr->pinning_ref_count <= 0 ) {

		unpin_entry(pinned_entry_ptr->cache_ptr,
			    pinned_entry_ptr->type,
			    pinned_entry_ptr->index);
	    }

	    entry_ptr->pin_type[i] = -1;
	    entry_ptr->pin_idx[i] = -1;
	}
	entry_ptr->num_pins = 0;
    }

    entry_ptr->destroyed = TRUE;
    entry_ptr->cache_ptr = NULL;

    return(SUCCEED);

} /* dest() */

herr_t
pico_dest(H5F_t * f, void * thing)
{
    HDassert ( ((test_entry_t *)thing)->type == PICO_ENTRY_TYPE );
    return(destroy(f, thing));
}

herr_t
nano_dest(H5F_t * f, void * thing)
{
    HDassert ( ((test_entry_t *)thing)->type == NANO_ENTRY_TYPE );
    return(destroy(f, thing));
}

herr_t
micro_dest(H5F_t * f, void * thing)
{
    HDassert ( ((test_entry_t *)thing)->type == MICRO_ENTRY_TYPE );
    return(destroy(f, thing));
}

herr_t
tiny_dest(H5F_t * f, void * thing)
{
    HDassert ( ((test_entry_t *)thing)->type == TINY_ENTRY_TYPE );
    return(destroy(f, thing));
}

herr_t
small_dest(H5F_t * f, void * thing)
{
    HDassert ( ((test_entry_t *)thing)->type == SMALL_ENTRY_TYPE );
    return(destroy(f, thing));
}

herr_t
medium_dest(H5F_t * f, void * thing)
{
    HDassert ( ((test_entry_t *)thing)->type == MEDIUM_ENTRY_TYPE );
    return(destroy(f, thing));
}

herr_t
large_dest(H5F_t * f, void * thing)
{
    HDassert ( ((test_entry_t *)thing)->type == LARGE_ENTRY_TYPE );
    return(destroy(f, thing));
}

herr_t
huge_dest(H5F_t * f, void *  thing)
{
    HDassert ( ((test_entry_t *)thing)->type == HUGE_ENTRY_TYPE );
    return(destroy(f, thing));
}

herr_t
monster_dest(H5F_t * f, void *  thing)
{
    HDassert ( ((test_entry_t *)thing)->type == MONSTER_ENTRY_TYPE );
    return(destroy(f, thing));
}

herr_t
variable_dest(H5F_t * f, void *  thing)
{
    HDassert ( ((test_entry_t *)thing)->type == VARIABLE_ENTRY_TYPE );
    return(destroy(f, thing));
}


/*-------------------------------------------------------------------------
 * Function:	flush & friends
 *
 * Purpose:	flush the entry and mark it as clean.  The helper functions
 *              verify that the correct version of flush is being called,
 *		and then call flush proper.
 *
 * Return:	SUCCEED
 *
 * Programmer:	John Mainzer
 *              6/10/04
 *
 * Modifications:
 *
 * 		JRM -- 8/30/06
 * 		Added variable_flush() and flags_ptr parameter.
 *
 * 		JRM -- 9/1/06
 * 		Added support for flush operations.
 *
 *-------------------------------------------------------------------------
 */

herr_t
flush(H5F_t *f,
      hid_t UNUSED dxpl_id,
      hbool_t dest,
      haddr_t
#ifdef NDEBUG
          UNUSED
#endif /* NDEBUG */
          addr,
      void *thing,
      unsigned * flags_ptr)
{
    int i;
    test_entry_t * entry_ptr;
    test_entry_t * base_addr;

    HDassert( thing );

    entry_ptr = (test_entry_t *)thing;
    base_addr = entries[entry_ptr->type];

    HDassert( entry_ptr->index >= 0 );
    HDassert( entry_ptr->index <= max_indices[entry_ptr->type] );
    HDassert( entry_ptr == &(base_addr[entry_ptr->index]) );
    HDassert( entry_ptr == entry_ptr->self );
    HDassert( entry_ptr->header.addr == entry_ptr->addr );
    HDassert( entry_ptr->addr == addr );
    HDassert( entry_ptr->header.size == entry_ptr->size );
    HDassert( ( entry_ptr->type == VARIABLE_ENTRY_TYPE ) ||
	      ( entry_ptr->size == entry_sizes[entry_ptr->type] ) );
    HDassert( entry_ptr->header.is_dirty == entry_ptr->is_dirty );
    HDassert( entry_ptr->cache_ptr != NULL );
    HDassert( entry_ptr->cache_ptr->magic == H5C__H5C_T_MAGIC );
    HDassert( entry_ptr->num_flush_ops >= 0 );
    HDassert( entry_ptr->num_flush_ops < MAX_FLUSH_OPS );

    if ( entry_ptr->num_flush_ops > 0 ) {

        for ( i = 0; i < entry_ptr->num_flush_ops; i++ )
	{
            execute_flush_op(entry_ptr->cache_ptr,
			     entry_ptr,
			     &((entry_ptr->flush_ops)[i]),
			     flags_ptr);
	}
	entry_ptr->num_flush_ops = 0;
	entry_ptr->flush_op_self_resize_in_progress = FALSE;
    }

    entry_ptr->flushed = TRUE;

    if ( ( ! write_permitted ) && ( entry_ptr->is_dirty ) ) {

        pass = FALSE;
        failure_mssg = "called flush when write_permitted is FALSE.";
    }

    if ( entry_ptr->is_dirty ) {

        (entry_ptr->writes)++;
        entry_ptr->is_dirty = FALSE;
        entry_ptr->header.is_dirty = FALSE;
    }

    if ( dest ) {

        destroy(f, thing);

    }

    return(SUCCEED);

} /* flush() */

herr_t
pico_flush(H5F_t *f, hid_t dxpl_id, hbool_t dest, haddr_t addr,
           void *thing, unsigned * flags_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == PICO_ENTRY_TYPE );
    return(flush(f, dxpl_id, dest, addr, thing, flags_ptr));
}

herr_t
nano_flush(H5F_t *f, hid_t dxpl_id, hbool_t dest, haddr_t addr,
	   void *thing, unsigned * flags_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == NANO_ENTRY_TYPE );
    return(flush(f, dxpl_id, dest, addr, thing, flags_ptr));
}

herr_t
micro_flush(H5F_t *f, hid_t dxpl_id, hbool_t dest, haddr_t addr,
            void *thing, unsigned * flags_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == MICRO_ENTRY_TYPE );
    return(flush(f, dxpl_id, dest, addr, thing, flags_ptr));
}

herr_t
tiny_flush(H5F_t *f, hid_t dxpl_id, hbool_t dest, haddr_t addr,
           void *thing, unsigned * flags_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == TINY_ENTRY_TYPE );
    return(flush(f, dxpl_id, dest, addr, thing, flags_ptr));
}

herr_t
small_flush(H5F_t *f, hid_t dxpl_id, hbool_t dest, haddr_t addr,
            void *thing, unsigned * flags_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == SMALL_ENTRY_TYPE );
    return(flush(f, dxpl_id, dest, addr, thing, flags_ptr));
}

herr_t
medium_flush(H5F_t *f, hid_t dxpl_id, hbool_t dest, haddr_t addr,
             void *thing, unsigned * flags_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == MEDIUM_ENTRY_TYPE );
    return(flush(f, dxpl_id, dest, addr, thing, flags_ptr));
}

herr_t
large_flush(H5F_t *f, hid_t dxpl_id, hbool_t dest, haddr_t addr,
            void *thing, unsigned * flags_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == LARGE_ENTRY_TYPE );
    return(flush(f, dxpl_id, dest, addr, thing, flags_ptr));
}

herr_t
huge_flush(H5F_t *f, hid_t dxpl_id, hbool_t dest, haddr_t addr,
           void *thing, unsigned * flags_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == HUGE_ENTRY_TYPE );
    return(flush(f, dxpl_id, dest, addr, thing, flags_ptr));
}

herr_t
monster_flush(H5F_t *f, hid_t dxpl_id, hbool_t dest, haddr_t addr,
	      void *thing, unsigned * flags_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == MONSTER_ENTRY_TYPE );
    return(flush(f, dxpl_id, dest, addr, thing, flags_ptr));
}

herr_t
variable_flush(H5F_t *f, hid_t dxpl_id, hbool_t dest, haddr_t addr,
	       void *thing, unsigned * flags_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == VARIABLE_ENTRY_TYPE );
    return(flush(f, dxpl_id, dest, addr, thing, flags_ptr));
}



/*-------------------------------------------------------------------------
 * Function:	load & friends
 *
 * Purpose:	"load" the requested entry and mark it as clean.  The
 *		helper functions verify that the correct version of load
 *		 is being called, and then call load proper.
 *
 * Return:	SUCCEED
 *
 * Programmer:	John Mainzer
 *              6/10/04
 *
 * Modifications:
 *
 * 		JRM -- 8/30/06
 * 		Added variable_load().
 *
 *-------------------------------------------------------------------------
 */

void *
load(H5F_t UNUSED *f,
     hid_t UNUSED dxpl_id,
     haddr_t addr,
     const void UNUSED *udata1,
     void UNUSED *udata2)
{
    int32_t type;
    int32_t idx;
    test_entry_t * entry_ptr;
    test_entry_t * base_addr;

    addr_to_type_and_index(addr, &type, &idx);

    base_addr = entries[type];
    entry_ptr = &(base_addr[idx]);

    HDassert( entry_ptr->type == type );
    HDassert( entry_ptr->type >= 0 );
    HDassert( entry_ptr->type < NUMBER_OF_ENTRY_TYPES );
    HDassert( entry_ptr->index == idx );
    HDassert( entry_ptr->index >= 0 );
    HDassert( entry_ptr->index <= max_indices[type] );
    HDassert( entry_ptr == entry_ptr->self );
    HDassert( entry_ptr->addr == addr );
#if 1 /* JRM */
    if ( ! ( ( entry_ptr->type == VARIABLE_ENTRY_TYPE ) ||
             ( entry_ptr->size == entry_sizes[type] ) ) ) {

        HDfprintf(stdout, "entry type/index/size = %d/%d/%ld\n",
                  (int)(entry_ptr->type),
                  (int)(entry_ptr->index),
                  (long)(entry_ptr->size));
    }
#endif /* JRM */
    HDassert( ( entry_ptr->type == VARIABLE_ENTRY_TYPE ) ||
	      ( entry_ptr->size == entry_sizes[type] ) );

    entry_ptr->loaded = TRUE;

    entry_ptr->header.is_dirty = FALSE;
    entry_ptr->is_dirty = FALSE;

    (entry_ptr->reads)++;

    return(entry_ptr);

} /* load() */

void *
pico_load(H5F_t *f, hid_t dxpl_id, haddr_t addr,
          const void *udata1, void *udata2)
{
    return(load(f, dxpl_id, addr, udata1, udata2));
}

void *
nano_load(H5F_t *f, hid_t dxpl_id, haddr_t addr,
          const void *udata1, void *udata2)
{
    return(load(f, dxpl_id, addr, udata1, udata2));
}

void *
micro_load(H5F_t *f, hid_t dxpl_id, haddr_t addr,
           const void *udata1, void *udata2)
{
    return(load(f, dxpl_id, addr, udata1, udata2));
}

void *
tiny_load(H5F_t *f, hid_t dxpl_id, haddr_t addr,
          const void *udata1, void *udata2)
{
    return(load(f, dxpl_id, addr, udata1, udata2));
}

void *
small_load(H5F_t *f, hid_t dxpl_id, haddr_t addr,
           const void *udata1, void *udata2)
{
    return(load(f, dxpl_id, addr, udata1, udata2));
}

void *
medium_load(H5F_t *f, hid_t dxpl_id, haddr_t addr,
            const void *udata1, void *udata2)
{
    return(load(f, dxpl_id, addr, udata1, udata2));
}

void *
large_load(H5F_t *f, hid_t dxpl_id, haddr_t addr,
           const void *udata1, void *udata2)
{
    return(load(f, dxpl_id, addr, udata1, udata2));
}

void *
huge_load(H5F_t *f, hid_t dxpl_id, haddr_t addr,
          const void *udata1, void *udata2)
{
    return(load(f, dxpl_id, addr, udata1, udata2));
}

void *
monster_load(H5F_t *f, hid_t dxpl_id, haddr_t addr,
             const void *udata1, void *udata2)
{
    return(load(f, dxpl_id, addr, udata1, udata2));
}

void *
variable_load(H5F_t *f, hid_t dxpl_id, haddr_t addr,
              const void *udata1, void *udata2)
{
    return(load(f, dxpl_id, addr, udata1, udata2));
}


/*-------------------------------------------------------------------------
 * Function:	size & friends
 *
 * Purpose:	Get the size of the specified entry.  The helper functions
 *		verify that the correct version of size is being called,
 *		and then call size proper.
 *
 * Return:	SUCCEED
 *
 * Programmer:	John Mainzer
 *              6/10/04
 *
 * Modifications:
 *
 * 		JRM -- 8/30/06
 * 		Added variable_size().
 *
 *-------------------------------------------------------------------------
 */

herr_t
size(H5F_t UNUSED *  f,
     void *   thing,
     size_t * size_ptr)
{
    test_entry_t * entry_ptr;
    test_entry_t * base_addr;

    HDassert( size_ptr );
    HDassert( thing );

    entry_ptr = (test_entry_t *)thing;
    base_addr = entries[entry_ptr->type];

    HDassert( entry_ptr->index >= 0 );
    HDassert( entry_ptr->index <= max_indices[entry_ptr->type] );
    HDassert( entry_ptr == &(base_addr[entry_ptr->index]) );
    HDassert( entry_ptr == entry_ptr->self );
    HDassert( entry_ptr->header.addr == entry_ptr->addr );
    HDassert( ( entry_ptr->type == VARIABLE_ENTRY_TYPE ) || \
              ( entry_ptr->size == entry_sizes[entry_ptr->type] ) );

    *size_ptr = entry_ptr->size;

    return(SUCCEED);

} /* size() */

herr_t
pico_size(H5F_t * f, void * thing, size_t * size_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == PICO_ENTRY_TYPE );
    return(size(f, thing, size_ptr));
}

herr_t
nano_size(H5F_t * f, void * thing, size_t * size_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == NANO_ENTRY_TYPE );
    return(size(f, thing, size_ptr));
}

herr_t
micro_size(H5F_t * f, void * thing, size_t * size_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == MICRO_ENTRY_TYPE );
    return(size(f, thing, size_ptr));
}

herr_t
tiny_size(H5F_t * f, void * thing, size_t * size_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == TINY_ENTRY_TYPE );
    return(size(f, thing, size_ptr));
}

herr_t
small_size(H5F_t * f, void * thing, size_t * size_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == SMALL_ENTRY_TYPE );
    return(size(f, thing, size_ptr));
}

herr_t
medium_size(H5F_t * f, void * thing, size_t * size_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == MEDIUM_ENTRY_TYPE );
    return(size(f, thing, size_ptr));
}

herr_t
large_size(H5F_t * f, void * thing, size_t * size_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == LARGE_ENTRY_TYPE );
    return(size(f, thing, size_ptr));
}

herr_t
huge_size(H5F_t * f, void * thing, size_t * size_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == HUGE_ENTRY_TYPE );
    return(size(f, thing, size_ptr));
}

herr_t
monster_size(H5F_t * f, void * thing, size_t * size_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == MONSTER_ENTRY_TYPE );
    return(size(f, thing, size_ptr));
}

herr_t
variable_size(H5F_t * f, void * thing, size_t * size_ptr)
{
    HDassert ( ((test_entry_t *)thing)->type == VARIABLE_ENTRY_TYPE );
    return(size(f, thing, size_ptr));
}



/**************************************************************************/
/**************************************************************************/
/************************** test utility functions: ***********************/
/**************************************************************************/
/**************************************************************************/

/*-------------------------------------------------------------------------
 * Function:	add_flush_op
 *
 * Purpose:	Do nothing if pass is FALSE on entry.
 *
 *              Otherwise, add the specified flush operation to the
 *              target instance of test_entry_t.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              9/1/06
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */

void
add_flush_op(int target_type,
	     int target_idx,
	     int op_code,
	     int type,
	     int idx,
	     hbool_t flag,
	     size_t new_size)
{
    int i;
    test_entry_t * target_base_addr;
    test_entry_t * target_entry_ptr;

    HDassert( ( 0 <= target_type ) && ( target_type < NUMBER_OF_ENTRY_TYPES ) );
    HDassert( ( 0 <= target_idx ) &&
	      ( target_idx <= max_indices[target_type] ) );
    HDassert( ( 0 <= op_code ) && ( op_code <= FLUSH_OP__MAX_OP ) );
    HDassert( ( op_code != FLUSH_OP__RESIZE ) ||
	      ( type == VARIABLE_ENTRY_TYPE ) );
    HDassert( ( 0 <= type ) && ( type < NUMBER_OF_ENTRY_TYPES ) );
    HDassert( ( 0 <= idx ) && ( idx <= max_indices[type] ) );
    HDassert( ( flag == TRUE ) || ( flag == FALSE ) );
    HDassert( new_size <= VARIABLE_ENTRY_SIZE );

    if ( pass ) {

        target_base_addr = entries[target_type];
        target_entry_ptr = &(target_base_addr[target_idx]);

        HDassert( target_entry_ptr->index == target_idx );
        HDassert( target_entry_ptr->type == target_type );
        HDassert( target_entry_ptr == target_entry_ptr->self );
	HDassert( target_entry_ptr->num_flush_ops < MAX_FLUSH_OPS );

	i = (target_entry_ptr->num_flush_ops)++;
	(target_entry_ptr->flush_ops)[i].op_code = op_code;
	(target_entry_ptr->flush_ops)[i].type = type;
	(target_entry_ptr->flush_ops)[i].idx = idx;
	(target_entry_ptr->flush_ops)[i].flag = flag;
	(target_entry_ptr->flush_ops)[i].size = new_size;

    }

    return;

} /* add_flush_op() */


/*-------------------------------------------------------------------------
 * Function:	create_pinned_entry_dependency
 *
 * Purpose:	Do nothing if pass is FALSE on entry.
 *
 *              Otherwise, set up a pinned entry dependency so we can
 *              test the pinned entry modifications to the flush routine.
 *
 *		Given the types and indicies of the pinned and pinning
 *		entries, add the pinned entry to the list of pinned
 *		entries in the pinning entry, increment the
 *		pinning reference count of the pinned entry, and
 *		if that count was zero initially, pin the entry.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              6/10/04
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */

void
create_pinned_entry_dependency(H5C_t * cache_ptr,
		               int pinning_type,
                               int pinning_idx,
	                       int pinned_type,
	                       int pinned_idx)
{
    test_entry_t * pinning_base_addr;
    test_entry_t * pinning_entry_ptr;
    test_entry_t * pinned_base_addr;
    test_entry_t * pinned_entry_ptr;

    if ( pass ) {

        HDassert( ( 0 <= pinning_type ) &&
 	          ( pinning_type < NUMBER_OF_ENTRY_TYPES ) );
        HDassert( ( 0 <= pinning_idx ) &&
	          ( pinning_idx <= max_indices[pinning_type] ) );
        HDassert( ( 0 <= pinned_type ) &&
	          ( pinned_type < NUMBER_OF_ENTRY_TYPES ) );
        HDassert( ( 0 <= pinned_idx ) &&
	          ( pinned_idx <= max_indices[pinned_type] ) );

        pinning_base_addr = entries[pinning_type];
        pinning_entry_ptr = &(pinning_base_addr[pinning_idx]);

        pinned_base_addr = entries[pinned_type];
        pinned_entry_ptr = &(pinned_base_addr[pinned_idx]);

        HDassert( pinning_entry_ptr->index == pinning_idx );
        HDassert( pinning_entry_ptr->type == pinning_type );
        HDassert( pinning_entry_ptr == pinning_entry_ptr->self );
	HDassert( pinning_entry_ptr->num_pins < MAX_PINS );

        HDassert( pinning_entry_ptr->index == pinning_idx );
        HDassert( pinning_entry_ptr->type == pinning_type );
        HDassert( pinning_entry_ptr == pinning_entry_ptr->self );
	HDassert( ! ( pinning_entry_ptr->is_protected ) );

	pinning_entry_ptr->pin_type[pinning_entry_ptr->num_pins] = pinned_type;
	pinning_entry_ptr->pin_idx[pinning_entry_ptr->num_pins] = pinned_idx;
	(pinning_entry_ptr->num_pins)++;

        if ( pinned_entry_ptr->pinning_ref_count == 0 ) {

	    protect_entry(cache_ptr, pinned_type, pinned_idx);
	    unprotect_entry(cache_ptr, pinned_type, pinned_idx, FALSE,
		            H5C__PIN_ENTRY_FLAG);
	}

	(pinned_entry_ptr->pinning_ref_count)++;
    }

    return;

} /* create_pinned_entry_dependency() */


/*-------------------------------------------------------------------------
 * Function:	dirty_entry
 *
 * Purpose:	Given a pointer to a cache, an entry type, and an index,
 *		dirty the target entry.
 *
 *		If the dirty_pin parameter is true, verify that the
 *		target entry is in the cache and is pinned.  If it
 *		isn't, scream and die.  If it is, use the
 *		H5C_mark_pinned_entry_dirty() call to dirty it.
 *
 *		Do nothing if pass is false on entry.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              6/10/04
 *
 * Modifications:
 *
 *		None.
 *
 *-------------------------------------------------------------------------
 */

void
dirty_entry(H5C_t * cache_ptr,
            int32_t type,
            int32_t idx,
	    hbool_t dirty_pin)
{
    test_entry_t * base_addr;
    test_entry_t * entry_ptr;

    HDassert( cache_ptr );
    HDassert( ( 0 <= type ) && ( type < NUMBER_OF_ENTRY_TYPES ) );
    HDassert( ( 0 <= idx ) && ( idx <= max_indices[type] ) );

    if ( pass ) {

        if ( dirty_pin ) {

	    if ( ! entry_in_cache(cache_ptr, type, idx) ) {

		pass = FALSE;
                failure_mssg = "entry to be dirty pinned is not in cache.";

	    } else {

                base_addr = entries[type];
                entry_ptr = &(base_addr[idx]);

	        HDassert( entry_ptr->index == idx );
	        HDassert( entry_ptr->type == type );
                HDassert( entry_ptr == entry_ptr->self );

		if ( ! ( (entry_ptr->header).is_pinned ) ) {

                    pass = FALSE;
                    failure_mssg = "entry to be dirty pinned is not pinned.";

                } else {

		    mark_pinned_entry_dirty(cache_ptr, type, idx, FALSE, (size_t)0);

		}
	    }
        } else {

	    protect_entry(cache_ptr, type, idx);
            unprotect_entry(cache_ptr, type, idx, TRUE, H5C__NO_FLAGS_SET);
	}
    }

    return;

} /* dirty_entry() */


/*-------------------------------------------------------------------------
 * Function:	execute_flush_op
 *
 * Purpose:	Given a pointer to an instance of struct flush_op, execute
 * 		it.
 *
 *		Do nothing if pass is false on entry.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              9/1/06
 *
 * Modifications:
 *
 *		None.
 *
 *-------------------------------------------------------------------------
 */

void
execute_flush_op(H5C_t * cache_ptr,
		 struct test_entry_t * entry_ptr,
		 struct flush_op * op_ptr,
		 unsigned * flags_ptr)
{
    HDassert( cache_ptr != NULL );
    HDassert( cache_ptr->magic == H5C__H5C_T_MAGIC );
    HDassert( entry_ptr != NULL );
    HDassert( entry_ptr = entry_ptr->self );
    HDassert( entry_ptr->header.addr == entry_ptr->addr );
    HDassert( ( entry_ptr->flush_op_self_resize_in_progress ) ||
              ( entry_ptr->header.size == entry_ptr->size ) );
    HDassert( op_ptr != NULL );
    HDassert( ( 0 <= entry_ptr->type ) &&
	      ( entry_ptr->type < NUMBER_OF_ENTRY_TYPES ) );
    HDassert( ( 0 <= entry_ptr->index ) &&
              ( entry_ptr->index <= max_indices[entry_ptr->type] ) );
    HDassert( ( 0 <= op_ptr->type ) &&
              ( op_ptr->type < NUMBER_OF_ENTRY_TYPES ) );
    HDassert( ( 0 <= op_ptr->idx ) &&
              ( op_ptr->idx <= max_indices[op_ptr->type] ) );
    HDassert( ( op_ptr->flag == FALSE ) || ( op_ptr->flag == TRUE ) );
    HDassert( flags_ptr != NULL );

    if ( pass ) {

	switch ( op_ptr->op_code )
	{
	    case FLUSH_OP__NO_OP:
		break;

	    case FLUSH_OP__DIRTY:
		HDassert( ( entry_ptr->type != op_ptr->type ) ||
			  ( entry_ptr->index != op_ptr->idx ) );

		dirty_entry(cache_ptr, op_ptr->type, op_ptr->idx, op_ptr->flag);
		break;

            case FLUSH_OP__RESIZE:
		if ( ( entry_ptr->type == op_ptr->type ) &&
                     ( entry_ptr->index == op_ptr->idx ) ) {

                    /* the flush operation is acting on the entry to
		     * which it is attached.  Handle this here:
		     */
                    HDassert( entry_ptr->type == VARIABLE_ENTRY_TYPE );
		    HDassert( op_ptr->size > 0 );
		    HDassert( op_ptr->size <= VARIABLE_ENTRY_SIZE );

                    entry_ptr->size = op_ptr->size;
		    (*flags_ptr) |= H5C_CALLBACK__SIZE_CHANGED_FLAG;
		    entry_ptr->flush_op_self_resize_in_progress = TRUE;

		    /* if the entry is in the process of being destroyed,
		     * set the header size to match the entry size so as
		     * to avoid a spurious failure in the destroy callback.
		     */
		    if ( entry_ptr->header.destroy_in_progress ) {

			entry_ptr->header.size = entry_ptr->size;
		    }

		} else {

		    /* change the size of some other entry */

		    resize_entry(cache_ptr, op_ptr->type, op_ptr->idx,
                                 op_ptr->size, op_ptr->flag);
		}
		break;

	    case FLUSH_OP__RENAME:
		rename_entry(cache_ptr, op_ptr->type, op_ptr->idx,
			     op_ptr->flag);
		break;

	    default:
                pass = FALSE;
                failure_mssg = "Undefined flush op code.";
		break;
	}
    }

    return;

} /* execute_flush_op() */


/*-------------------------------------------------------------------------
 * Function:	entry_in_cache
 *
 * Purpose:	Given a pointer to a cache, an entry type, and an index,
 *		determine if the entry is currently in the cache.
 *
 * Return:	TRUE if the entry is in the cache, and FALSE otherwise.
 *
 * Programmer:	John Mainzer
 *              6/10/04
 *
 * Modifications:
 *
 *		JRM - 10/12/04
 *		Removed references to local_H5C_t, as we now get direct
 *		access to the definition of H5C_t via H5Cpkg.h.
 *
 *-------------------------------------------------------------------------
 */

hbool_t
entry_in_cache(H5C_t * cache_ptr,
               int32_t type,
               int32_t idx)
{
    hbool_t in_cache = FALSE; /* will set to TRUE if necessary */
    test_entry_t * base_addr;
    test_entry_t * entry_ptr;
    H5C_cache_entry_t * test_ptr = NULL;

    HDassert( cache_ptr );
    HDassert( ( 0 <= type ) && ( type < NUMBER_OF_ENTRY_TYPES ) );
    HDassert( ( 0 <= idx ) && ( idx <= max_indices[type] ) );

    base_addr = entries[type];
    entry_ptr = &(base_addr[idx]);

    HDassert( entry_ptr->index == idx );
    HDassert( entry_ptr->type == type );
    HDassert( entry_ptr == entry_ptr->self );

    H5C__SEARCH_INDEX(cache_ptr, entry_ptr->addr, test_ptr)

    if ( test_ptr != NULL ) {

        in_cache = TRUE;
        HDassert( test_ptr == (H5C_cache_entry_t *)entry_ptr );
        HDassert( entry_ptr->addr == entry_ptr->header.addr );
    }

    return(in_cache);

} /* entry_in_cache() */


/*-------------------------------------------------------------------------
 * Function:	reset_entries
 *
 * Purpose:	reset the contents of the entries arrays to know values.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              6/10/04
 *
 * Modifications:
 *
 * 		JRM -- 3/31/06
 * 		Added initialization for new pinned entry test related
 * 		fields.
 *
 * 		JRM -- 4/1/07
 * 		Added initialization for the new is_read_only, and
 * 		ro_ref_count fields.
 *
 *-------------------------------------------------------------------------
 */

void
reset_entries(void)

{
    int i;
    int j;
    int k;
    int32_t max_index;
    haddr_t addr = 0;
    haddr_t alt_addr = PICO_ALT_BASE_ADDR;
    size_t entry_size;
    test_entry_t * base_addr;

    for ( i = 0; i < NUMBER_OF_ENTRY_TYPES; i++ )
    {
        entry_size = entry_sizes[i];
        max_index = max_indices[i];
        base_addr = entries[i];

        HDassert( base_addr );

        for ( j = 0; j <= max_index; j++ )
        {
            /* one can argue that we should fill the header with garbage.
             * If this is desired, we can simply comment out the header
             * initialization - the headers will be full of garbage soon
             * enough.
             */

            base_addr[j].header.addr = (haddr_t)0;
            base_addr[j].header.size = (size_t)0;
            base_addr[j].header.type = NULL;
            base_addr[j].header.is_dirty = FALSE;
            base_addr[j].header.is_protected = FALSE;
            base_addr[j].header.is_read_only = FALSE;
            base_addr[j].header.ro_ref_count = FALSE;
            base_addr[j].header.next = NULL;
            base_addr[j].header.prev = NULL;
            base_addr[j].header.aux_next = NULL;
            base_addr[j].header.aux_prev = NULL;

            base_addr[j].self = &(base_addr[j]);
            base_addr[j].cache_ptr = NULL;
            base_addr[j].addr = addr;
            base_addr[j].at_main_addr = TRUE;
            base_addr[j].main_addr = addr;
            base_addr[j].alt_addr = alt_addr;
            base_addr[j].size = entry_size;
            base_addr[j].type = i;
            base_addr[j].index = j;
            base_addr[j].reads = 0;
            base_addr[j].writes = 0;
            base_addr[j].is_dirty = FALSE;
            base_addr[j].is_protected = FALSE;
            base_addr[j].is_read_only = FALSE;
            base_addr[j].ro_ref_count = FALSE;

            base_addr[j].is_pinned = FALSE;
	    base_addr[j].pinning_ref_count = 0;
	    base_addr[j].num_pins = 0;
	    for ( k = 0; k < MAX_PINS; k++ )
            {
	        base_addr[j].pin_type[k] = -1;
		base_addr[j].pin_idx[k] = -1;
	    }

	    base_addr[j].num_flush_ops = 0;
	    for ( k = 0; k < MAX_FLUSH_OPS; k++ )
	    {
		base_addr[j].flush_ops[k].op_code = FLUSH_OP__NO_OP;
		base_addr[j].flush_ops[k].type = -1;
		base_addr[j].flush_ops[k].idx = -1;
		base_addr[j].flush_ops[k].flag = FALSE;
		base_addr[j].flush_ops[k].size = 0;
            }
	    base_addr[j].flush_op_self_resize_in_progress = FALSE;

            base_addr[j].loaded = FALSE;
            base_addr[j].cleared = FALSE;
            base_addr[j].flushed = FALSE;
            base_addr[j].destroyed = FALSE;

            addr += (haddr_t)entry_size;
            alt_addr += (haddr_t)entry_size;
        }
    }

    return;

} /* reset_entries() */


/*-------------------------------------------------------------------------
 * Function:	resize_entry
 *
 * Purpose:	Given a pointer to a cache, an entry type, an index, and
 * 		a size, set the size of the target entry to the size.  Note
 * 		that at present, the type of the entry must be
 * 		VARIABLE_ENTRY_TYPE.
 *
 *		If the resize_pin parameter is true, verify that the
 *		target entry is in the cache and is pinned.  If it
 *		isn't, scream and die.  If it is, use the
 *		H5C_mark_pinned_entry_dirty() call to resize it.
 *
 *		Do nothing if pass is false on entry.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              6/10/04
 *
 * Modifications:
 *
 *		None.
 *
 *-------------------------------------------------------------------------
 */

void
resize_entry(H5C_t * cache_ptr,
             int32_t type,
             int32_t idx,
	     size_t new_size,
	     hbool_t resize_pin)
{
    test_entry_t * base_addr;
    test_entry_t * entry_ptr;

    HDassert( cache_ptr );
    HDassert( ( 0 <= type ) && ( type < NUMBER_OF_ENTRY_TYPES ) );
    HDassert( type == VARIABLE_ENTRY_TYPE );
    HDassert( ( 0 <= idx ) && ( idx <= max_indices[type] ) );
    HDassert( ( 0 < new_size ) && ( new_size <= entry_sizes[type] ) );

    if ( pass ) {

        base_addr = entries[type];
        entry_ptr = &(base_addr[idx]);

        HDassert( entry_ptr->index == idx );
        HDassert( entry_ptr->type == type );
        HDassert( entry_ptr == entry_ptr->self );

        if ( resize_pin ) {

	    if ( ! entry_in_cache(cache_ptr, type, idx) ) {

		pass = FALSE;
                failure_mssg = "entry to be resized pinned is not in cache.";

	    } else {

		if ( ! ( (entry_ptr->header).is_pinned ) ) {

                    pass = FALSE;
                    failure_mssg = "entry to be resized pinned is not pinned.";

                } else {

		    mark_pinned_entry_dirty(cache_ptr, type, idx,
				            TRUE, new_size);
		}
	    }
        } else {

	    protect_entry(cache_ptr, type, idx);
	    unprotect_entry_with_size_change(cache_ptr, type, idx,
                                             H5C__SIZE_CHANGED_FLAG, new_size);
	}
    }

    return;

} /* resize_entry() */


/*-------------------------------------------------------------------------
 * Function:	resize_pinned_entry
 *
 * Purpose:	Given a pointer to a cache, an entry type, an index, and
 *              a new size, change the size of the target pinned entry
 *              to match the supplied new size.
 *
 *		Do nothing if pass is false on entry.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              1/11/08
 *
 * Modifications:
 *
 *		None.
 *
 *-------------------------------------------------------------------------
 */

void
resize_pinned_entry(H5C_t * cache_ptr,
                    int32_t type,
                    int32_t idx,
	            size_t new_size)
{
    herr_t result;
    test_entry_t * base_addr;
    test_entry_t * entry_ptr;

    HDassert( cache_ptr );
    HDassert( ( 0 <= type ) && ( type < NUMBER_OF_ENTRY_TYPES ) );
    HDassert( ( 0 <= idx ) && ( idx <= max_indices[type] ) );
    HDassert( type == VARIABLE_ENTRY_TYPE ) ;
    HDassert( ( 0 < new_size ) && ( new_size <= entry_sizes[type] ) );

    if ( pass ) {

        if ( ! entry_in_cache(cache_ptr, type, idx) ) {

	    pass = FALSE;
            failure_mssg = "entry not in cache.";

        } else {

            base_addr = entries[type];
            entry_ptr = &(base_addr[idx]);

            HDassert( entry_ptr->index == idx );
            HDassert( entry_ptr->type == type );
            HDassert( entry_ptr == entry_ptr->self );

            if ( ! ( (entry_ptr->header).is_pinned ) ) {

                pass = FALSE;
                failure_mssg = "entry to be resized is not pinned.";

            } else {

		entry_ptr->size = new_size;

	        result = H5C_resize_pinned_entry(cache_ptr,
				                 (void *)entry_ptr,
						 new_size);

		if ( result != SUCCEED ) {

		    pass = FALSE;
		    failure_mssg = "error(s) in H5C_resize_pinned_entry().";

		} else {

		    HDassert( entry_ptr->size = (entry_ptr->header).size );

                }
	    }
	}
    }

    return;

} /* resize_pinned_entry() */


/*-------------------------------------------------------------------------
 * Function:	verify_clean
 *
 * Purpose:	Verify that all cache entries are marked as clean.  If any
 *		are not, set pass to FALSE.
 *
 *		Do nothing if pass is FALSE on entry.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              6/10/04
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */

void
verify_clean(void)

{
    int i;
    int j;
    int dirty_count = 0;
    int32_t max_index;
    test_entry_t * base_addr;

    if ( pass ) {

        for ( i = 0; i < NUMBER_OF_ENTRY_TYPES; i++ )
        {
            max_index = max_indices[i];
            base_addr = entries[i];

            HDassert( base_addr );

            for ( j = 0; j <= max_index; j++ )
            {
                if ( ( base_addr[j].header.is_dirty ) ||
		     ( base_addr[j].is_dirty ) ) {

                    dirty_count++;
                }
            }
        }

        if ( dirty_count > 0 ) {

            pass = FALSE;
            failure_mssg = "verify_clean() found dirty entry(s).";
        }
    }

    return;

} /* verify_clean() */


/*-------------------------------------------------------------------------
 * Function:	verify_entry_status
 *
 * Purpose:	Verify that a list of entries have the expected status.
 * 		If any discrepencies are found, set the failure message
 * 		and set pass to FALSE.
 *
 *		Do nothing if pass is FALSE on entry.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              10/8/04
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */

void
verify_entry_status(H5C_t * cache_ptr,
		    int tag,
		    int num_entries,
		    struct expected_entry_status expected[])
{
    static char    msg[128];
    hbool_t        in_cache = FALSE; /* will set to TRUE if necessary */
    int            i;
    test_entry_t * entry_ptr;
    test_entry_t * base_addr;

    i = 0;
    while ( ( pass ) && ( i < num_entries ) )
    {
        base_addr = entries[expected[i].entry_type];
	entry_ptr = &(base_addr[expected[i].entry_index]);

	if ( ( ! expected[i].in_cache ) &&
	     ( ( expected[i].is_dirty ) ||
	       ( expected[i].is_protected ) ||
	       ( expected[i].is_pinned ) ) ) {

	    pass = FALSE;
	    sprintf(msg, "Contradictory data in expected[%d].\n", i);
	    failure_mssg = msg;
	}

        if ( pass ) {

	    in_cache = entry_in_cache(cache_ptr, expected[i].entry_type,
		                      expected[i].entry_index);

	    if ( in_cache != expected[i].in_cache ) {

	        pass = FALSE;
	        sprintf(msg,
		      "%d entry (%d, %d) in cache actual/expected = %d/%d.\n",
		      tag,
		      (int)expected[i].entry_type,
		      (int)expected[i].entry_index,
		      (int)in_cache,
		      (int)expected[i].in_cache);
	        failure_mssg = msg;
	    }
	}

        if ( pass ) {

	    if ( entry_ptr->size != expected[i].size ) {

	        pass = FALSE;
	        sprintf(msg,
                        "%d entry (%d, %d) size actualexpected = %ld/%ld.\n",
			tag,
	                (int)expected[i].entry_type,
		        (int)expected[i].entry_index,
		        (long)(entry_ptr->size),
		        (long)expected[i].size);
	        failure_mssg = msg;
	    }
	}

        if ( ( pass ) && ( in_cache ) ) {

	    if ( entry_ptr->header.size != expected[i].size ) {

	        pass = FALSE;
	        sprintf(msg,
                        "%d entry (%d, %d) header size actual/expected = %ld/%ld.\n",
			tag,
		        (int)expected[i].entry_type,
		        (int)expected[i].entry_index,
		        (long)(entry_ptr->header.size),
		        (long)expected[i].size);
	        failure_mssg = msg;
	    }
	}

	if ( pass ) {

	    if ( entry_ptr->at_main_addr != expected[i].at_main_addr ) {

	        pass = FALSE;
	        sprintf(msg,
                      "%d entry (%d, %d) at main addr actual/expected = %d/%d.\n",
		      tag,
		      (int)expected[i].entry_type,
		      (int)expected[i].entry_index,
		      (int)(entry_ptr->at_main_addr),
		      (int)expected[i].at_main_addr);
	        failure_mssg = msg;
	    }
	}

	if ( pass ) {

	    if ( entry_ptr->is_dirty != expected[i].is_dirty ) {

	        pass = FALSE;
	        sprintf(msg,
                      "%d entry (%d, %d) is_dirty actual/expected = %d/%d.\n",
		      tag,
		      (int)expected[i].entry_type,
		      (int)expected[i].entry_index,
		      (int)(entry_ptr->is_dirty),
		      (int)expected[i].is_dirty);
	        failure_mssg = msg;
	    }
	}

	if ( ( pass ) && ( in_cache ) ) {

	    if ( entry_ptr->header.is_dirty != expected[i].is_dirty ) {

	        pass = FALSE;
	        sprintf(msg,
                      "%d entry (%d, %d) header is_dirty actual/expected = %d/%d.\n",
		      tag,
		      (int)expected[i].entry_type,
		      (int)expected[i].entry_index,
		      (int)(entry_ptr->header.is_dirty),
		      (int)expected[i].is_dirty);
	        failure_mssg = msg;
	    }
	}

	if ( pass ) {

	    if ( entry_ptr->is_protected != expected[i].is_protected ) {

	        pass = FALSE;
	        sprintf(msg,
                      "%d entry (%d, %d) is_protected actual/expected = %d/%d.\n",
		      tag,
		      (int)expected[i].entry_type,
		      (int)expected[i].entry_index,
		      (int)(entry_ptr->is_protected),
		      (int)expected[i].is_protected);
	        failure_mssg = msg;
	    }
	}

	if ( ( pass ) && ( in_cache ) ) {

	    if ( entry_ptr->header.is_protected != expected[i].is_protected ) {

	        pass = FALSE;
	        sprintf(msg,
                      "%d entry (%d, %d) header is_protected actual/expected = %d/%d.\n",
		      tag,
		      (int)expected[i].entry_type,
		      (int)expected[i].entry_index,
		      (int)(entry_ptr->header.is_protected),
		      (int)expected[i].is_protected);
	        failure_mssg = msg;
	    }
	}

	if ( pass ) {

	    if ( entry_ptr->is_pinned != expected[i].is_pinned ) {

	        pass = FALSE;
	        sprintf(msg,
                      "%d entry (%d, %d) is_pinned actual/expected = %d/%d.\n",
		      tag,
		      (int)expected[i].entry_type,
		      (int)expected[i].entry_index,
		      (int)(entry_ptr->is_pinned),
		      (int)expected[i].is_pinned);
	        failure_mssg = msg;
	    }
	}

	if ( ( pass ) && ( in_cache ) ) {

	    if ( entry_ptr->header.is_pinned != expected[i].is_pinned ) {

	        pass = FALSE;
	        sprintf(msg,
                  "%d entry (%d, %d) header is_pinned actual/expected = %d/%d.\n",
		  tag,
		  (int)expected[i].entry_type,
		  (int)expected[i].entry_index,
		  (int)(entry_ptr->header.is_pinned),
		  (int)expected[i].is_pinned);
	        failure_mssg = msg;
	    }
	}

	if ( pass ) {

            if ( ( entry_ptr->loaded != expected[i].loaded ) ||
	         ( entry_ptr->cleared != expected[i].cleared ) ||
	         ( entry_ptr->flushed != expected[i].flushed ) ||
	         ( entry_ptr->destroyed != expected[i].destroyed ) ) {

	        pass = FALSE;
                sprintf(msg,
                        "%d entry (%d,%d) loaded = %d(%d), clrd = %d(%d), flshd = %d(%d), dest = %d(%d)\n",
			tag,
		        (int)expected[i].entry_type,
		        (int)expected[i].entry_index,
		        (int)(entry_ptr->loaded),
		        (int)(expected[i].loaded),
		        (int)(entry_ptr->cleared),
		        (int)(expected[i].cleared),
		        (int)(entry_ptr->flushed),
		        (int)(expected[i].flushed),
		        (int)(entry_ptr->destroyed),
		        (int)(expected[i].destroyed));
                failure_mssg = msg;
            }
        }
	i++;
    } /* while */

    return;

} /* verify_entry_status() */


/*-------------------------------------------------------------------------
 * Function:	verify_unprotected
 *
 * Purpose:	Verify that no cache entries are marked as protected.  If
 *		any are, set pass to FALSE.
 *
 *		Do nothing if pass is FALSE on entry.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              6/10/04
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */

void
verify_unprotected(void)

{
    int i;
    int j;
    int protected_count = 0;
    int32_t max_index;
    test_entry_t * base_addr;

    if ( pass ) {

        for ( i = 0; i < NUMBER_OF_ENTRY_TYPES; i++ )
        {
            max_index = max_indices[i];
            base_addr = entries[i];

            HDassert( base_addr );

            for ( j = 0; j <= max_index; j++ )
            {
                HDassert( base_addr[j].header.is_protected ==
                          base_addr[j].is_protected );

                if ( ( base_addr[j].header.is_protected ) ||
                     ( base_addr[j].is_protected ) ) {

                    protected_count++;
                }
            }
        }

        if ( protected_count > 0 ) {

            pass = FALSE;
            failure_mssg = "verify_unprotected() found protected entry(s).";
        }
    }

    return;

} /* verify_unprotected() */


/*-------------------------------------------------------------------------
 * Function:	setup_cache()
 *
 * Purpose:	Allocate a cache of the desired size and configure it for
 *		use in the test bed.  Return a pointer to the new cache
 *		structure.
 *
 * Return:	Pointer to new cache, or NULL on failure.
 *
 * Programmer:	John Mainzer
 *              6/11/04
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */

H5C_t *
setup_cache(size_t max_cache_size,
            size_t min_clean_size)
{
    H5C_t * cache_ptr = NULL;

    cache_ptr = H5C_create(max_cache_size,
                           min_clean_size,
                           (NUMBER_OF_ENTRY_TYPES - 1),
			   (const char **)entry_type_names,
                           check_write_permitted,
                           TRUE,
                           NULL,
                           NULL);

    if ( cache_ptr == NULL ) {

        pass = FALSE;
        failure_mssg = "H5C_create() returned NULL.";

    } else {

        H5C_set_skip_flags(cache_ptr, TRUE, TRUE);
    }

    return(cache_ptr);

} /* setup_cache() */


/*-------------------------------------------------------------------------
 * Function:	takedown_cache()
 *
 * Purpose:	Flush the specified cache and disable it.  If requested,
 *		dump stats first.  If pass is FALSE, do nothing.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              6/11/04
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */

void
takedown_cache(H5C_t * cache_ptr,
               hbool_t dump_stats,
               hbool_t dump_detailed_stats)
{
    HDassert(cache_ptr);

    if ( pass ) {

        if ( dump_stats ) {

            H5C_stats(cache_ptr, "test cache", dump_detailed_stats);
        }

        H5C_dest(NULL, -1, -1, cache_ptr);
    }

    return;

} /* takedown_cache() */


/*-------------------------------------------------------------------------
 * Function:	expunge_entry()
 *
 * Purpose:	Expunge the entry indicated by the type and index.
 *
 *		Do nothing if pass is FALSE on entry.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              7/6/06
 *
 * Modifications:
 *
 *		None.
 *
 *-------------------------------------------------------------------------
 */

void
expunge_entry(H5C_t * cache_ptr,
              int32_t type,
              int32_t idx)
{
    /* const char * fcn_name = "expunge_entry()"; */
    herr_t result;
    test_entry_t * base_addr;
    test_entry_t * entry_ptr;

    if ( pass ) {

        HDassert( cache_ptr );
        HDassert( ( 0 <= type ) && ( type < NUMBER_OF_ENTRY_TYPES ) );
        HDassert( ( 0 <= idx ) && ( idx <= max_indices[type] ) );

        base_addr = entries[type];
        entry_ptr = &(base_addr[idx]);

        HDassert( entry_ptr->index == idx );
        HDassert( entry_ptr->type == type );
        HDassert( entry_ptr == entry_ptr->self );
	HDassert( entry_ptr->cache_ptr == cache_ptr );
        HDassert( ! ( entry_ptr->header.is_protected ) );
        HDassert( ! ( entry_ptr->is_protected ) );
        HDassert( ! ( entry_ptr->header.is_pinned ) );
	HDassert( ! ( entry_ptr->is_pinned ) );

        result = H5C_expunge_entry(NULL, -1, -1, cache_ptr, &(types[type]),
                                   entry_ptr->addr, H5AC__NO_FLAGS_SET);

        if ( result < 0 ) {

            pass = FALSE;
            failure_mssg = "error in H5C_expunge_entry().";

        }
    }

    return;

} /* expunge_entry() */


/*-------------------------------------------------------------------------
 * Function:	flush_cache()
 *
 * Purpose:	Flush the specified cache, destroying all entries if
                requested.  If requested, dump stats first.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              6/23/04
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */

void
flush_cache(H5C_t * cache_ptr,
            hbool_t destroy_entries,
            hbool_t dump_stats,
            hbool_t dump_detailed_stats)
{
    const char * fcn_name = "flush_cache()";
    herr_t result = 0;
    hbool_t verbose = FALSE;

    HDassert(cache_ptr);

    verify_unprotected();

    if ( pass ) {

        if ( destroy_entries ) {

            result = H5C_flush_cache(NULL, -1, -1, cache_ptr,
                                     H5C__FLUSH_INVALIDATE_FLAG);

        } else {

            result = H5C_flush_cache(NULL, -1, -1, cache_ptr,
                                     H5C__NO_FLAGS_SET);
        }
    }

    if ( dump_stats ) {

        H5C_stats(cache_ptr, "test cache", dump_detailed_stats);
    }

    if ( result < 0 ) {

        pass = FALSE;
        failure_mssg = "error in H5C_flush_cache().";
    }
    else if ( ( destroy_entries ) &&
	      ( ( cache_ptr->index_len != 0 ) ||
	        ( cache_ptr->index_size != 0 ) ||
	        ( cache_ptr->clean_index_size != 0 ) ||
	        ( cache_ptr->dirty_index_size != 0 ) ) ) {

	if ( verbose ) {
            HDfprintf(stdout, 
		      "%s: unexpected il/is/cis/dis = %lld/%lld/%lld/%lld.\n",
		      fcn_name,
		      (long long)(cache_ptr->index_len),
		      (long long)(cache_ptr->index_size),
		      (long long)(cache_ptr->clean_index_size),
		      (long long)(cache_ptr->dirty_index_size));
	}
        pass = FALSE;
        failure_mssg = 
	   "non zero index len/sizes after H5C_flush_cache() with invalidate.";
    }


    return;

} /* flush_cache() */


/*-------------------------------------------------------------------------
 * Function:	insert_entry()
 *
 * Purpose:	Insert the entry indicated by the type and index.  Mark
 *		it clean or dirty as indicated.
 *
 *		Note that I don't see much practical use for inserting
 *		a clean entry, but the interface permits it so we should
 *		test it.
 *
 *		Do nothing if pass is false.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              6/16/04
 *
 * Modifications:
 *
 *		JRM -- 1/13/05
 *		Updated function for the flags parameter in
 *		H5C_insert_entry(), and to allow access to this parameter.
 *
 *		JRM -- 6/17/05
 *		The interface no longer permits clean inserts.
 *		Accordingly, the dirty parameter is no longer meaningfull.
 *
 *		JRM -- 4/5/06
 *		Added code to initialize the new cache_ptr field of the
 *		test_entry_t structure.
 *
 *		JRM -- 8/10/06
 *		Updated to reflect the fact that entries can now be
 *		inserted pinned.
 *
 *-------------------------------------------------------------------------
 */

void
insert_entry(H5C_t * cache_ptr,
             int32_t type,
             int32_t idx,
             hbool_t UNUSED dirty,
             unsigned int flags)
{
    herr_t result;
    hbool_t insert_pinned;
    test_entry_t * base_addr;
    test_entry_t * entry_ptr;

    if ( pass ) {

        HDassert( cache_ptr );
        HDassert( ( 0 <= type ) && ( type < NUMBER_OF_ENTRY_TYPES ) );
        HDassert( ( 0 <= idx ) && ( idx <= max_indices[type] ) );

        base_addr = entries[type];
        entry_ptr = &(base_addr[idx]);

        HDassert( entry_ptr->index == idx );
        HDassert( entry_ptr->type == type );
        HDassert( entry_ptr == entry_ptr->self );
        HDassert( !(entry_ptr->is_protected) );

	insert_pinned = ((flags & H5C__PIN_ENTRY_FLAG) != 0 );

	entry_ptr->is_dirty = TRUE;

        result = H5C_insert_entry(NULL, -1, -1, cache_ptr, &(types[type]),
                                  entry_ptr->addr, (void *)entry_ptr, flags);

        if ( ( result < 0 ) ||
             ( entry_ptr->header.is_protected ) ||
             ( entry_ptr->header.type != &(types[type]) ) ||
             ( entry_ptr->size != entry_ptr->header.size ) ||
             ( entry_ptr->addr != entry_ptr->header.addr ) ) {

            pass = FALSE;
            failure_mssg = "error in H5C_insert().";

#if 0 /* This is useful debugging code.  Lets keep it around. */

            HDfprintf(stdout, "result = %d\n", (int)result);
            HDfprintf(stdout, "entry_ptr->header.is_protected = %d\n",
                      (int)(entry_ptr->header.is_protected));
            HDfprintf(stdout,
		      "entry_ptr->header.type != &(types[type]) = %d\n",
                      (int)(entry_ptr->header.type != &(types[type])));
            HDfprintf(stdout,
                      "entry_ptr->size != entry_ptr->header.size = %d\n",
                      (int)(entry_ptr->size != entry_ptr->header.size));
            HDfprintf(stdout,
                      "entry_ptr->addr != entry_ptr->header.addr = %d\n",
                       (int)(entry_ptr->addr != entry_ptr->header.addr));
#endif
        }
	HDassert( entry_ptr->cache_ptr == NULL );

        entry_ptr->cache_ptr = cache_ptr;

	if ( insert_pinned ) {

	    HDassert( entry_ptr->header.is_pinned );
	    entry_ptr->is_pinned = TRUE;

	} else {

	    HDassert( ! ( entry_ptr->header.is_pinned ) );
	    entry_ptr->is_pinned = FALSE;

	}
        HDassert( entry_ptr->header.is_dirty );
        HDassert( ((entry_ptr->header).type)->id == type );
    }

    return;

} /* insert_entry() */


/*-------------------------------------------------------------------------
 * Function:	mark_pinned_entry_dirty()
 *
 * Purpose:	Mark the specified entry as dirty.
 *
 *		Do nothing if pass is FALSE on entry.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              3/28/06
 *
 * Modifications:
 *
 *		None.
 *
 *-------------------------------------------------------------------------
 */

void
mark_pinned_entry_dirty(H5C_t * cache_ptr,
                        int32_t type,
                        int32_t idx,
			hbool_t size_changed,
			size_t  new_size)
{
    /* const char * fcn_name = "mark_pinned_entry_dirty()"; */
    herr_t result;
    test_entry_t * base_addr;
    test_entry_t * entry_ptr;

    if ( pass ) {

        HDassert( cache_ptr );
        HDassert( ( 0 <= type ) && ( type < NUMBER_OF_ENTRY_TYPES ) );
        HDassert( ( 0 <= idx ) && ( idx <= max_indices[type] ) );

        base_addr = entries[type];
        entry_ptr = &(base_addr[idx]);

        HDassert( entry_ptr->index == idx );
        HDassert( entry_ptr->type == type );
        HDassert( entry_ptr == entry_ptr->self );
	HDassert( entry_ptr->cache_ptr == cache_ptr );
        HDassert( ! (entry_ptr->header.is_protected) );
        HDassert( entry_ptr->header.is_pinned );
	HDassert( entry_ptr->is_pinned );

	entry_ptr->is_dirty = TRUE;

	if ( size_changed ) {

	    /* update entry size now to keep the sanity checks happy */
	    entry_ptr->size = new_size;
	}

        result = H5C_mark_pinned_entry_dirty(cache_ptr,
			                     (void *)entry_ptr,
					     size_changed,
					     new_size);

        if ( ( result < 0 ) ||
             ( ! (entry_ptr->header.is_dirty) ) ||
             ( ! (entry_ptr->header.is_pinned) ) ||
             ( entry_ptr->header.type != &(types[type]) ) ||
             ( entry_ptr->size != entry_ptr->header.size ) ||
             ( entry_ptr->addr != entry_ptr->header.addr ) ) {

#if 0 /* This is useful debugging code -- keep it around  */
	    HDfprintf(stdout, "result = %ld.\n", (long)result);
	    HDfprintf(stdout, "entry_ptr->header.is_dirty = %d.\n",
		      (int)(entry_ptr->header.is_dirty));
	    HDfprintf(stdout, "entry_ptr->header.is_pinned = %d.\n",
		      (int)(entry_ptr->header.is_pinned));
	    HDfprintf(stdout,
		      "(entry_ptr->header.type != &(types[type])) = %d.\n",
		      (int)(entry_ptr->header.type != &(types[type])));
	    HDfprintf(stdout,
		      "entry_ptr->size = %ld, entry_ptr->header.size = %ld.\n",
		      (long)(entry_ptr->size), (long)(entry_ptr->header.size));
	    HDfprintf(stdout,
		      "entry_ptr->addr = %ld, entry_ptr->header.addr = %ld.\n",
		      (long)(entry_ptr->addr), (long)(entry_ptr->header.addr));
#endif
            pass = FALSE;
            failure_mssg = "error in H5C_mark_pinned_entry_dirty().";

        }

        HDassert( ((entry_ptr->header).type)->id == type );

    }

    return;

} /* mark_pinned_entry_dirty() */


/*-------------------------------------------------------------------------
 * Function:	mark_pinned_or_protected_entry_dirty()
 *
 * Purpose:	Mark the specified entry as dirty.
 *
 *		Do nothing if pass is FALSE on entry.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              5/17/06
 *
 * Modifications:
 *
 *		None.
 *
 *-------------------------------------------------------------------------
 */

void
mark_pinned_or_protected_entry_dirty(H5C_t * cache_ptr,
                                     int32_t type,
                                     int32_t idx)
{
    /* const char * fcn_name = "mark_pinned_or_protected_entry_dirty()"; */
    herr_t result;
    test_entry_t * base_addr;
    test_entry_t * entry_ptr;

    if ( pass ) {

        HDassert( cache_ptr );
        HDassert( ( 0 <= type ) && ( type < NUMBER_OF_ENTRY_TYPES ) );
        HDassert( ( 0 <= idx ) && ( idx <= max_indices[type] ) );

        base_addr = entries[type];
        entry_ptr = &(base_addr[idx]);

        HDassert( entry_ptr->index == idx );
        HDassert( entry_ptr->type == type );
        HDassert( entry_ptr == entry_ptr->self );
	HDassert( entry_ptr->cache_ptr == cache_ptr );
        HDassert( entry_ptr->header.is_protected ||
		  entry_ptr->header.is_pinned );

	entry_ptr->is_dirty = TRUE;

        result = H5C_mark_pinned_or_protected_entry_dirty(cache_ptr,
			                                  (void *)entry_ptr);

        if ( ( result < 0 )
	     ||
	     ( ( ! (entry_ptr->header.is_protected) )
	       &&
	       ( ! (entry_ptr->header.is_pinned) )
	     )
	     ||
             ( ( entry_ptr->header.is_protected )
	       &&
	       ( ! ( entry_ptr->header.dirtied ) )
	     )
	     ||
             ( ( ! ( entry_ptr->header.is_protected ) )
	       &&
	       ( ! ( entry_ptr->header.is_dirty ) )
	     )
	     ||
             ( entry_ptr->header.type != &(types[type]) )
	     ||
             ( entry_ptr->size != entry_ptr->header.size )
	     ||
             ( entry_ptr->addr != entry_ptr->header.addr ) ) {

            pass = FALSE;
            failure_mssg =
                "error in H5C_mark_pinned_or_protected_entry_dirty().";

        }

        HDassert( ((entry_ptr->header).type)->id == type );

    }

    return;

} /* mark_pinned_or_protected_entry_dirty() */


/*-------------------------------------------------------------------------
 * Function:	rename_entry()
 *
 * Purpose:	Rename the entry indicated by the type and index to its
 *		main or alternate address as indicated.  If the entry is
 *		already at the desired entry, do nothing.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              6/21/04
 *
 * Modifications:
 *
 *		JRM -- 6/17/05
 *		Updated code to reflect the fact that renames automatically
 *		dirty entries.
 *
 *-------------------------------------------------------------------------
 */

void
rename_entry(H5C_t * cache_ptr,
             int32_t type,
             int32_t idx,
             hbool_t main_addr)
{
    herr_t         result;
    hbool_t	   done = TRUE; /* will set to FALSE if we have work to do */
    haddr_t        old_addr = HADDR_UNDEF;
    haddr_t        new_addr = HADDR_UNDEF;
    test_entry_t * base_addr;
    test_entry_t * entry_ptr;

    HDassert( cache_ptr );
    HDassert( ( 0 <= type ) && ( type < NUMBER_OF_ENTRY_TYPES ) );
    HDassert( ( 0 <= idx ) && ( idx <= max_indices[type] ) );

    base_addr = entries[type];
    entry_ptr = &(base_addr[idx]);

    HDassert( entry_ptr->index == idx );
    HDassert( entry_ptr->type == type );
    HDassert( entry_ptr == entry_ptr->self );
    HDassert( entry_ptr->cache_ptr == cache_ptr );
    HDassert( !(entry_ptr->is_protected) );
    HDassert( !(entry_ptr->header.is_protected) );


    if ( entry_ptr->at_main_addr && !main_addr ) {

        /* rename to alt addr */

        HDassert( entry_ptr->addr == entry_ptr->main_addr );

        done = FALSE;
        old_addr = entry_ptr->addr;
        new_addr = entry_ptr->alt_addr;

    } else if ( !(entry_ptr->at_main_addr) && main_addr ) {

        /* rename to main addr */

        HDassert( entry_ptr->addr == entry_ptr->alt_addr );

        done = FALSE;
        old_addr = entry_ptr->addr;
        new_addr = entry_ptr->main_addr;
    }

    if ( ! done ) {

        entry_ptr->is_dirty = TRUE;

        result = H5C_rename_entry(cache_ptr, &(types[type]),
                                  old_addr, new_addr);
    }

    if ( ! done ) {

        if ( ( result < 0 ) ||
	     ( ( ! ( entry_ptr->header.destroy_in_progress ) ) &&
	       ( entry_ptr->header.addr != new_addr ) ) ) {

            pass = FALSE;
            failure_mssg = "error in H5C_rename_entry().";

        } else {

            entry_ptr->addr = new_addr;
            entry_ptr->at_main_addr = main_addr;
        }
    }

    HDassert( ((entry_ptr->header).type)->id == type );

    HDassert( entry_ptr->header.is_dirty );
    HDassert( entry_ptr->is_dirty );

    return;

} /* rename_entry() */


/*-------------------------------------------------------------------------
 * Function:	protect_entry()
 *
 * Purpose:	Protect the entry indicated by the type and index.
 *
 *		Do nothing if pass is FALSE on entry.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              6/11/04
 *
 * Modifications:
 *
 *    - Modified call to H5C_protect to pass H5C__NO_FLAGS_SET in the
 *      new flags parameter.
 *    						JRM -- 3/28/07
 *
 *-------------------------------------------------------------------------
 */

void
protect_entry(H5C_t * cache_ptr,
              int32_t type,
              int32_t idx)
{
    /* const char * fcn_name = "protect_entry()"; */
    test_entry_t * base_addr;
    test_entry_t * entry_ptr;
    H5C_cache_entry_t * cache_entry_ptr;

    if ( pass ) {

        HDassert( cache_ptr );
        HDassert( ( 0 <= type ) && ( type < NUMBER_OF_ENTRY_TYPES ) );
        HDassert( ( 0 <= idx ) && ( idx <= max_indices[type] ) );

        base_addr = entries[type];
        entry_ptr = &(base_addr[idx]);

        HDassert( entry_ptr->index == idx );
        HDassert( entry_ptr->type == type );
        HDassert( entry_ptr == entry_ptr->self );
        HDassert( !(entry_ptr->is_protected) );

        cache_entry_ptr = H5C_protect(NULL, -1, -1, cache_ptr, &(types[type]),
                                      entry_ptr->addr, NULL, NULL,
				      H5C__NO_FLAGS_SET);

        if ( ( cache_entry_ptr != (void *)entry_ptr ) ||
             ( !(entry_ptr->header.is_protected) ) ||
             ( entry_ptr->header.type != &(types[type]) ) ||
             ( entry_ptr->size != entry_ptr->header.size ) ||
             ( entry_ptr->addr != entry_ptr->header.addr ) ) {

#if 0
            /* I've written the following debugging code several times
             * now.  Lets keep it around so I don't have to write it
             * again.
             *                              - JRM
             */
            HDfprintf(stdout, "( cache_entry_ptr != (void *)entry_ptr ) = %d\n",
                      (int)( cache_entry_ptr != (void *)entry_ptr ));
            HDfprintf(stdout, "cache_entry_ptr = 0x%lx, entry_ptr = 0x%lx\n",
                      (long)cache_entry_ptr, (long)entry_ptr);
            HDfprintf(stdout, "entry_ptr->header.is_protected = %d\n",
                      (int)(entry_ptr->header.is_protected));
            HDfprintf(stdout,
                      "( entry_ptr->header.type != &(types[type]) ) = %d\n",
                      (int)( entry_ptr->header.type != &(types[type]) ));
            HDfprintf(stdout,
                      "entry_ptr->size = %d, entry_ptr->header.size = %d\n",
                      (int)(entry_ptr->size), (int)(entry_ptr->header.size));
            HDfprintf(stdout,
                      "entry_ptr->addr = %d, entry_ptr->header.addr = %d\n",
                      (int)(entry_ptr->addr), (int)(entry_ptr->header.addr));
#endif
            pass = FALSE;
            failure_mssg = "error in H5C_protect().";

        } else {

	    HDassert( ( entry_ptr->cache_ptr == NULL ) ||
		      ( entry_ptr->cache_ptr == cache_ptr ) );

	    entry_ptr->cache_ptr = cache_ptr;
            entry_ptr->is_protected = TRUE;

        }

        HDassert( ((entry_ptr->header).type)->id == type );
    }

    return;

} /* protect_entry() */


/*-------------------------------------------------------------------------
 * Function:	protect_entry_ro()
 *
 * Purpose:	Do a read only protect the entry indicated by the type
 * 		and index.
 *
 *		Do nothing if pass is FALSE on entry.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              4/1/07
 *
 * Modifications:
 *
 *    - None.
 *
 *-------------------------------------------------------------------------
 */

void
protect_entry_ro(H5C_t * cache_ptr,
                int32_t type,
                int32_t idx)
{
    /* const char * fcn_name = "protect_entry_ro()"; */
    test_entry_t * base_addr;
    test_entry_t * entry_ptr;
    H5C_cache_entry_t * cache_entry_ptr;

    if ( pass ) {

        HDassert( cache_ptr );
        HDassert( ( 0 <= type ) && ( type < NUMBER_OF_ENTRY_TYPES ) );
        HDassert( ( 0 <= idx ) && ( idx <= max_indices[type] ) );

        base_addr = entries[type];
        entry_ptr = &(base_addr[idx]);

        HDassert( entry_ptr->index == idx );
        HDassert( entry_ptr->type == type );
        HDassert( entry_ptr == entry_ptr->self );
        HDassert( ( ! ( entry_ptr->is_protected ) ) ||
		  ( ( entry_ptr->is_read_only ) &&
		    ( entry_ptr->ro_ref_count > 0 ) ) );

        cache_entry_ptr = H5C_protect(NULL, -1, -1, cache_ptr, &(types[type]),
                                      entry_ptr->addr, NULL, NULL,
				      H5C__READ_ONLY_FLAG);

        if ( ( cache_entry_ptr != (void *)entry_ptr ) ||
             ( !(entry_ptr->header.is_protected) ) ||
             ( !(entry_ptr->header.is_read_only) ) ||
             ( entry_ptr->header.ro_ref_count <= 0 ) ||
             ( entry_ptr->header.type != &(types[type]) ) ||
             ( entry_ptr->size != entry_ptr->header.size ) ||
             ( entry_ptr->addr != entry_ptr->header.addr ) ) {

            pass = FALSE;
            failure_mssg = "error in read only H5C_protect().";

        } else {

	    HDassert( ( entry_ptr->cache_ptr == NULL ) ||
		      ( entry_ptr->cache_ptr == cache_ptr ) );

	    entry_ptr->cache_ptr = cache_ptr;
            entry_ptr->is_protected = TRUE;
	    entry_ptr->is_read_only = TRUE;
	    entry_ptr->ro_ref_count++;
        }

        HDassert( ((entry_ptr->header).type)->id == type );
    }

    return;

} /* protect_entry_ro() */


/*-------------------------------------------------------------------------
 * Function:	unpin_entry()
 *
 * Purpose:	Unpin the entry indicated by the type and index.
 *
 *		Do nothing if pass is FALSE on entry.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              3/28/06
 *
 * Modifications:
 *
 *		None.
 *
 *-------------------------------------------------------------------------
 */

void
unpin_entry(H5C_t * cache_ptr,
            int32_t type,
            int32_t idx)
{
    /* const char * fcn_name = "unpin_entry()"; */
    herr_t result;
    test_entry_t * base_addr;
    test_entry_t * entry_ptr;

    if ( pass ) {

        HDassert( cache_ptr );
        HDassert( ( 0 <= type ) && ( type < NUMBER_OF_ENTRY_TYPES ) );
        HDassert( ( 0 <= idx ) && ( idx <= max_indices[type] ) );

        base_addr = entries[type];
        entry_ptr = &(base_addr[idx]);

        HDassert( entry_ptr->index == idx );
        HDassert( entry_ptr->type == type );
        HDassert( entry_ptr == entry_ptr->self );
	HDassert( entry_ptr->cache_ptr == cache_ptr );
        HDassert( ! (entry_ptr->header.is_protected) );
        HDassert( entry_ptr->header.is_pinned );
	HDassert( entry_ptr->is_pinned );

        result = H5C_unpin_entry(cache_ptr, (void *)entry_ptr);

        if ( ( result < 0 ) ||
             ( entry_ptr->header.is_pinned ) ||
             ( entry_ptr->header.type != &(types[type]) ) ||
             ( entry_ptr->size != entry_ptr->header.size ) ||
             ( entry_ptr->addr != entry_ptr->header.addr ) ) {

            pass = FALSE;
            failure_mssg = "error in H5C_unpin().";

        }

	entry_ptr->is_pinned = FALSE;

        HDassert( ((entry_ptr->header).type)->id == type );

    }

    return;

} /* unpin_entry() */


/*-------------------------------------------------------------------------
 * Function:	unprotect_entry()
 *
 * Purpose:	Unprotect the entry indicated by the type and index.
 *
 *		Do nothing if pass is FALSE on entry.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              6/12/04
 *
 * Modifications:
 *
 *		JRM -- 1/7/05
 *		Updated for the replacement of the deleted parameter in
 *		H5C_unprotect() with the new flags parameter.
 *
 *		JRM - 6/17/05
 *		Modified function to use the new dirtied parameter of
 *		H5C_unprotect().
 *
 *		JRM -- 9/8/05
 *		Update for new entry size parameter in H5C_unprotect().
 *		We don't use them here for now.
 *
 *		JRM -- 3/31/06
 *		Update for pinned entries.
 *
 *		JRM -- 4/1/07
 *		Updated for new multiple read protects.
 *
 *-------------------------------------------------------------------------
 */

void
unprotect_entry(H5C_t * cache_ptr,
                int32_t type,
                int32_t idx,
                int dirty,
                unsigned int flags)
{
    /* const char * fcn_name = "unprotect_entry()"; */
    herr_t result;
    hbool_t pin_flag_set;
    hbool_t unpin_flag_set;
    test_entry_t * base_addr;
    test_entry_t * entry_ptr;

    if ( pass ) {

        HDassert( cache_ptr );
        HDassert( ( 0 <= type ) && ( type < NUMBER_OF_ENTRY_TYPES ) );
        HDassert( ( 0 <= idx ) && ( idx <= max_indices[type] ) );

        base_addr = entries[type];
        entry_ptr = &(base_addr[idx]);

        HDassert( entry_ptr->index == idx );
        HDassert( entry_ptr->type == type );
        HDassert( entry_ptr == entry_ptr->self );
	HDassert( entry_ptr->cache_ptr == cache_ptr );
        HDassert( entry_ptr->header.is_protected );
        HDassert( entry_ptr->is_protected );

	pin_flag_set = ((flags & H5C__PIN_ENTRY_FLAG) != 0 );
	unpin_flag_set = ((flags & H5C__UNPIN_ENTRY_FLAG) != 0 );

	HDassert ( ! ( pin_flag_set && unpin_flag_set ) );
	HDassert ( ( ! pin_flag_set ) || ( ! (entry_ptr->is_pinned) ) );
	HDassert ( ( ! unpin_flag_set ) || ( entry_ptr->is_pinned ) );

        if ( ( dirty == TRUE ) || ( dirty == FALSE ) ) {

            flags |= (dirty ? H5C__DIRTIED_FLAG : H5C__NO_FLAGS_SET);
            entry_ptr->is_dirty = (entry_ptr->is_dirty || dirty);
        }

        result = H5C_unprotect(NULL, -1, -1, cache_ptr, &(types[type]),
                               entry_ptr->addr, (void *)entry_ptr,
                               flags, (size_t)0);

        if ( ( result < 0 ) ||
             ( ( entry_ptr->header.is_protected ) &&
	       ( ( ! ( entry_ptr->is_read_only ) ) ||
		 ( entry_ptr->ro_ref_count <= 0 ) ) ) ||
             ( entry_ptr->header.type != &(types[type]) ) ||
             ( entry_ptr->size != entry_ptr->header.size ) ||
             ( entry_ptr->addr != entry_ptr->header.addr ) ) {

#if 1 /* JRM */
	    if ( result < 0 ) {
		HDfprintf(stdout, "result is negative.\n");
	    }
	    if ( ( entry_ptr->header.is_protected ) &&
                 ( ( ! ( entry_ptr->is_read_only ) ) ||
                   ( entry_ptr->ro_ref_count <= 0 ) ) ) {
		HDfprintf(stdout, "protected and not RO or refcnt <= 0.\n");
	    }
            if ( entry_ptr->header.type != &(types[type]) ) {
		HDfprintf(stdout, "type disagreement.\n");
	    }
	    if ( entry_ptr->size != entry_ptr->header.size ) {
		HDfprintf(stdout, "size disagreement.\n");
	    }
	    if ( entry_ptr->addr != entry_ptr->header.addr ) {
		HDfprintf(stdout, "addr disagreement.\n");
	    }
#endif /* JRM */

            pass = FALSE;
            failure_mssg = "error in H5C_unprotect().";

        }
        else
        {
	    if ( entry_ptr->ro_ref_count > 1 ) {

		entry_ptr->ro_ref_count--;

	    } else if ( entry_ptr->ro_ref_count == 1 ) {

		entry_ptr->is_protected = FALSE;
		entry_ptr->is_read_only = FALSE;
		entry_ptr->ro_ref_count = 0;

	    } else {

		entry_ptr->is_protected = FALSE;

	    }

	    if ( pin_flag_set ) {

	        HDassert ( entry_ptr->header.is_pinned );
		entry_ptr->is_pinned = TRUE;

	    } else if ( unpin_flag_set ) {

	        HDassert ( ! ( entry_ptr->header.is_pinned ) );
		entry_ptr->is_pinned = FALSE;

            }
        }

        HDassert( ((entry_ptr->header).type)->id == type );

        if ( ( flags & H5C__DIRTIED_FLAG ) != 0
                && ( (flags & H5C__DELETED_FLAG) == 0 ) ) {

            HDassert( entry_ptr->header.is_dirty );
            HDassert( entry_ptr->is_dirty );
        }

	HDassert( entry_ptr->header.is_protected == entry_ptr->is_protected );
	HDassert( entry_ptr->header.is_read_only == entry_ptr->is_read_only );
	HDassert( entry_ptr->header.ro_ref_count == entry_ptr->ro_ref_count );
    }

    return;

} /* unprotect_entry() */


/*-------------------------------------------------------------------------
 * Function:	unprotect_entry_with_size_change()
 *
 * Purpose:	Version of unprotect_entry() that allow access to the new
 * 		size change parameters in H5C_unprotect_entry()
 *
 * 		At present, only the sizes of VARIABLE_ENTRY_TYPE entries
 * 		can be changed.  Thus this function will scream and die
 * 		if the H5C__SIZE_CHANGED_FLAG is set and the type is not
 * 		VARIABLE_ENTRY_TYPE.
 *
 *		Do nothing if pass is FALSE on entry.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              8/31/06
 *
 * Modifications:
 *
 *		None.
 *
 *-------------------------------------------------------------------------
 */

void
unprotect_entry_with_size_change(H5C_t * cache_ptr,
                                 int32_t type,
                                 int32_t idx,
                                 unsigned int flags,
		                 size_t new_size)
{
    const char * fcn_name = "unprotect_entry_with_size_change()";
    herr_t result;
    hbool_t dirty_flag_set;
    hbool_t pin_flag_set;
    hbool_t unpin_flag_set;
    hbool_t size_changed_flag_set;
    hbool_t verbose = FALSE;
    test_entry_t * base_addr;
    test_entry_t * entry_ptr;

    if ( pass ) {

        HDassert( cache_ptr );
        HDassert( ( 0 <= type ) && ( type < NUMBER_OF_ENTRY_TYPES ) );
        HDassert( ( 0 <= idx ) && ( idx <= max_indices[type] ) );
	HDassert( new_size <= entry_sizes[type] );

        base_addr = entries[type];
        entry_ptr = &(base_addr[idx]);

        HDassert( entry_ptr->index == idx );
        HDassert( entry_ptr->type == type );
        HDassert( entry_ptr == entry_ptr->self );
	HDassert( entry_ptr->cache_ptr == cache_ptr );
        HDassert( entry_ptr->header.is_protected );
        HDassert( entry_ptr->is_protected );

	dirty_flag_set = ((flags & H5C__DIRTIED_FLAG) != 0 );
	pin_flag_set = ((flags & H5C__PIN_ENTRY_FLAG) != 0 );
	unpin_flag_set = ((flags & H5C__UNPIN_ENTRY_FLAG) != 0 );
	size_changed_flag_set = ((flags & H5C__SIZE_CHANGED_FLAG) != 0 );

	HDassert ( ! ( pin_flag_set && unpin_flag_set ) );
	HDassert ( ( ! pin_flag_set ) || ( ! (entry_ptr->is_pinned) ) );
	HDassert ( ( ! unpin_flag_set ) || ( entry_ptr->is_pinned ) );
	HDassert ( ( ! size_changed_flag_set ) || ( new_size > 0 ) );
	HDassert ( ( ! size_changed_flag_set ) ||
		   ( type == VARIABLE_ENTRY_TYPE ) );

        entry_ptr->is_dirty = (entry_ptr->is_dirty || dirty_flag_set);

	if ( size_changed_flag_set ) {

            entry_ptr->is_dirty = TRUE;
	    entry_ptr->size = new_size;
        }

        result = H5C_unprotect(NULL, -1, -1, cache_ptr, &(types[type]),
                               entry_ptr->addr, (void *)entry_ptr,
                               flags, new_size);

        if ( ( result < 0 ) ||
             ( entry_ptr->header.is_protected ) ||
             ( entry_ptr->header.type != &(types[type]) ) ||
             ( entry_ptr->size != entry_ptr->header.size ) ||
             ( entry_ptr->addr != entry_ptr->header.addr ) ) {

	    if ( verbose ) {

	        if ( result < 0 ) {
                    HDfprintf(stdout, "%s: H5C_unprotect() failed.\n", fcn_name);
		}

	        if ( entry_ptr->header.is_protected )  {
                    HDfprintf(stdout, "%s: entry still protected?!?.\n", 
			      fcn_name);
		}

	        if ( entry_ptr->header.type != &(types[type]) )  {
                    HDfprintf(stdout, 
			      "%s: entry has bad type after unprotect.\n", 
			      fcn_name);
		}

	        if ( entry_ptr->size != entry_ptr->header.size )  {
                    HDfprintf(stdout, 
			   "%s: bad entry size after unprotect. e/a = %d/%d\n", 
			   fcn_name,
			   (int)(entry_ptr->size),
			   (int)(entry_ptr->header.size));
		}

	        if ( entry_ptr->addr != entry_ptr->header.addr )  {
                    HDfprintf(stdout, 
		     "%s: bad entry addr after unprotect. e/a = 0x%llx/0x%llx\n",
		     fcn_name,
		     (long long)(entry_ptr->addr),
		     (long long)(entry_ptr->header.addr));
		}
	    }
            pass = FALSE;
            failure_mssg = "error in H5C_unprotect().";

        }
        else
        {
            entry_ptr->is_protected = FALSE;

	    if ( pin_flag_set ) {

	        HDassert ( entry_ptr->header.is_pinned );
		entry_ptr->is_pinned = TRUE;

	    } else if ( unpin_flag_set ) {

	        HDassert ( ! ( entry_ptr->header.is_pinned ) );
		entry_ptr->is_pinned = FALSE;

            }
        }

        HDassert( ((entry_ptr->header).type)->id == type );

        if ( ( flags & H5C__DIRTIED_FLAG ) != 0
                && ( (flags & H5C__DELETED_FLAG) == 0 ) ) {

            HDassert( entry_ptr->header.is_dirty );
            HDassert( entry_ptr->is_dirty );
        }
    }

    return;

} /* unprotect_entry_with_size_change() */


/*-------------------------------------------------------------------------
 * Function:	row_major_scan_forward()
 *
 * Purpose:	Do a sequence of inserts, protects, unprotects, renames,
 *		destroys while scanning through the set of entries.  If
 *		pass is false on entry, do nothing.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              6/12/04
 *
 * Modifications:
 *
 * 		JRM -- 4/4/07
 * 		Added code supporting multiple read only protects.
 * 		Note that this increased the minimum lag to 10.
 *
 *-------------------------------------------------------------------------
 */

void
row_major_scan_forward(H5C_t * cache_ptr,
                       int32_t lag,
                       hbool_t verbose,
                       hbool_t reset_stats,
                       hbool_t display_stats,
                       hbool_t display_detailed_stats,
                       hbool_t do_inserts,
                       hbool_t dirty_inserts,
                       hbool_t do_renames,
                       hbool_t rename_to_main_addr,
                       hbool_t do_destroys,
		       hbool_t do_mult_ro_protects,
                       int dirty_destroys,
                       int dirty_unprotects)
{
    const char * fcn_name = "row_major_scan_forward";
    int32_t type;
    int32_t idx;

    if ( verbose )
        HDfprintf(stdout, "%s(): entering.\n", fcn_name);

    HDassert( lag >= 10 );

    type = 0;

    if ( ( pass ) && ( reset_stats ) ) {

        H5C_stats__reset(cache_ptr);
    }

    while ( ( pass ) && ( type < NUMBER_OF_ENTRY_TYPES ) )
    {
        idx = -lag;

        while ( ( pass ) && ( idx <= (max_indices[type] + lag) ) )
        {
	    if ( verbose ) {

                HDfprintf(stdout, "%d:%d: ", type, idx);
	    }

            if ( ( pass ) && ( do_inserts ) && ( (idx + lag) >= 0 ) &&
                 ( (idx + lag) <= max_indices[type] ) &&
                 ( ((idx + lag) % 2) == 0 ) &&
                 ( ! entry_in_cache(cache_ptr, type, (idx + lag)) ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(i, %d, %d) ", type, (idx + lag));

                insert_entry(cache_ptr, type, (idx + lag), dirty_inserts,
                             H5C__NO_FLAGS_SET);
            }


            if ( ( pass ) && ( (idx + lag - 1) >= 0 ) &&
                 ( (idx + lag - 1) <= max_indices[type] ) &&
                 ( ( (idx + lag - 1) % 3 ) == 0 ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(p, %d, %d) ", type, (idx + lag - 1));

                protect_entry(cache_ptr, type, (idx + lag - 1));
            }

            if ( ( pass ) && ( (idx + lag - 2) >= 0 ) &&
                 ( (idx + lag - 2) <= max_indices[type] ) &&
                 ( ( (idx + lag - 2) % 3 ) == 0 ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(u, %d, %d) ", type, (idx + lag - 2));

                unprotect_entry(cache_ptr, type, idx+lag-2, NO_CHANGE,
                                H5C__NO_FLAGS_SET);
            }


            if ( ( pass ) && ( do_renames ) && ( (idx + lag - 2) >= 0 ) &&
                 ( (idx + lag - 2) <= max_indices[type] ) &&
                 ( ( (idx + lag - 2) % 3 ) == 0 ) ) {

                rename_entry(cache_ptr, type, (idx + lag - 2),
                             rename_to_main_addr);
            }


            if ( ( pass ) && ( (idx + lag - 3) >= 0 ) &&
                 ( (idx + lag - 3) <= max_indices[type] ) &&
                 ( ( (idx + lag - 3) % 5 ) == 0 ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(p, %d, %d) ", type, (idx + lag - 3));

                protect_entry(cache_ptr, type, (idx + lag - 3));
            }

            if ( ( pass ) && ( (idx + lag - 5) >= 0 ) &&
                 ( (idx + lag - 5) <= max_indices[type] ) &&
                 ( ( (idx + lag - 5) % 5 ) == 0 ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(u, %d, %d) ", type, (idx + lag - 5));

                unprotect_entry(cache_ptr, type, idx+lag-5, NO_CHANGE,
                                H5C__NO_FLAGS_SET);
            }

	    if ( do_mult_ro_protects )
	    {
		if ( ( pass ) && ( (idx + lag - 5) >= 0 ) &&
		     ( (idx + lag - 5) < max_indices[type] ) &&
		     ( (idx + lag - 5) % 9 == 0 ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(p-ro, %d, %d) ", type,
				  (idx + lag - 5));

		    protect_entry_ro(cache_ptr, type, (idx + lag - 5));
		}

		if ( ( pass ) && ( (idx + lag - 6) >= 0 ) &&
		     ( (idx + lag - 6) < max_indices[type] ) &&
		     ( (idx + lag - 6) % 11 == 0 ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(p-ro, %d, %d) ", type,
				  (idx + lag - 6));

		    protect_entry_ro(cache_ptr, type, (idx + lag - 6));
		}

		if ( ( pass ) && ( (idx + lag - 7) >= 0 ) &&
		     ( (idx + lag - 7) < max_indices[type] ) &&
		     ( (idx + lag - 7) % 13 == 0 ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(p-ro, %d, %d) ", type,
				  (idx + lag - 7));

		    protect_entry_ro(cache_ptr, type, (idx + lag - 7));
		}

		if ( ( pass ) && ( (idx + lag - 7) >= 0 ) &&
		     ( (idx + lag - 7) < max_indices[type] ) &&
		     ( (idx + lag - 7) % 9 == 0 ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(u-ro, %d, %d) ", type,
				  (idx + lag - 7));

		    unprotect_entry(cache_ptr, type, (idx + lag - 7),
				    FALSE, H5C__NO_FLAGS_SET);
		}

		if ( ( pass ) && ( (idx + lag - 8) >= 0 ) &&
		     ( (idx + lag - 8) < max_indices[type] ) &&
		     ( (idx + lag - 8) % 11 == 0 ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(u-ro, %d, %d) ", type,
				  (idx + lag - 8));

		    unprotect_entry(cache_ptr, type, (idx + lag - 8),
				    FALSE, H5C__NO_FLAGS_SET);
		}

		if ( ( pass ) && ( (idx + lag - 9) >= 0 ) &&
		     ( (idx + lag - 9) < max_indices[type] ) &&
		     ( (idx + lag - 9) % 13 == 0 ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(u-ro, %d, %d) ", type,
				  (idx + lag - 9));

		    unprotect_entry(cache_ptr, type, (idx + lag - 9),
				    FALSE, H5C__NO_FLAGS_SET);
		}
	    } /* if ( do_mult_ro_protects ) */

            if ( ( pass ) && ( idx >= 0 ) && ( idx <= max_indices[type] ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(p, %d, %d) ", type, idx);

                protect_entry(cache_ptr, type, idx);
            }

            if ( ( pass ) && ( (idx - lag + 2) >= 0 ) &&
                 ( (idx - lag + 2) <= max_indices[type] ) &&
                 ( ( (idx - lag + 2) % 7 ) == 0 ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(u, %d, %d) ", type, (idx - lag + 2));

                unprotect_entry(cache_ptr, type, idx-lag+2, NO_CHANGE,
                                H5C__NO_FLAGS_SET);
            }

            if ( ( pass ) && ( (idx - lag + 1) >= 0 ) &&
                 ( (idx - lag + 1) <= max_indices[type] ) &&
                 ( ( (idx - lag + 1) % 7 ) == 0 ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(p, %d, %d) ", type, (idx - lag + 1));

                protect_entry(cache_ptr, type, (idx - lag + 1));
            }


            if ( do_destroys ) {

                if ( ( pass ) && ( (idx - lag) >= 0 ) &&
                     ( ( idx - lag) <= max_indices[type] ) ) {

                    switch ( (idx - lag) %4 ) {

                        case 0: /* we just did an insert */
                            unprotect_entry(cache_ptr, type, idx - lag,
                                            NO_CHANGE, H5C__NO_FLAGS_SET);
                            break;

                        case 1:
                            if ( (entries[type])[idx-lag].is_dirty ) {

                                unprotect_entry(cache_ptr, type, idx - lag,
                                                NO_CHANGE, H5C__NO_FLAGS_SET);
                            } else {

                                unprotect_entry(cache_ptr, type, idx - lag,
                                                dirty_unprotects,
                                                H5C__NO_FLAGS_SET);
                            }
                            break;

                        case 2: /* we just did an insrt */
                            unprotect_entry(cache_ptr, type, idx - lag,
                                            NO_CHANGE, H5C__DELETED_FLAG);
                            break;

                        case 3:
                            if ( (entries[type])[idx-lag].is_dirty ) {

                                unprotect_entry(cache_ptr, type, idx - lag,
                                                NO_CHANGE, H5C__DELETED_FLAG);
                            } else {

                                unprotect_entry(cache_ptr, type, idx - lag,
                                                dirty_destroys,
                                                H5C__DELETED_FLAG);
                            }
                            break;

                        default:
                            HDassert(0); /* this can't happen... */
                            break;
                    }
                }

            } else {

                if ( ( pass ) && ( (idx - lag) >= 0 ) &&
                     ( ( idx - lag) <= max_indices[type] ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(u, %d, %d) ", type, (idx - lag));

                    unprotect_entry(cache_ptr, type, idx - lag,
                                    dirty_unprotects, H5C__NO_FLAGS_SET);
                }
            }

            if ( verbose )
                HDfprintf(stdout, "\n");

            idx++;
        }
        type++;
    }

    if ( ( pass ) && ( display_stats ) ) {

        H5C_stats(cache_ptr, "test cache", display_detailed_stats);
    }

    return;

} /* row_major_scan_forward() */


/*-------------------------------------------------------------------------
 * Function:	hl_row_major_scan_forward()
 *
 * Purpose:	Do a high locality sequence of inserts, protects, and
 *		unprotects while scanning through the set of entries.
 *		If pass is false on entry, do nothing.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              10/21/04
 *
 * Modifications:
 *
 *		JRM -- 1/21/05
 *		Added the max_index parameter to allow the caller to
 *		throttle the size of the inner loop, and thereby the
 *		execution time of the function.
 *
 *-------------------------------------------------------------------------
 */

void
hl_row_major_scan_forward(H5C_t * cache_ptr,
                          int32_t max_index,
                          hbool_t verbose,
                          hbool_t reset_stats,
                          hbool_t display_stats,
                          hbool_t display_detailed_stats,
                          hbool_t do_inserts,
                          hbool_t dirty_inserts)
{
    const char * fcn_name = "hl_row_major_scan_forward";
    int32_t type;
    int32_t idx;
    int32_t i;
    int32_t lag = 100;
    int32_t local_max_index;

    if ( verbose )
        HDfprintf(stdout, "%s(): entering.\n", fcn_name);

    HDassert( lag > 5 );
    HDassert( max_index >= 200 );
    HDassert( max_index <= MAX_ENTRIES );

    type = 0;

    if ( ( pass ) && ( reset_stats ) ) {

        H5C_stats__reset(cache_ptr);
    }

    while ( ( pass ) && ( type < NUMBER_OF_ENTRY_TYPES ) )
    {
        idx = -lag;

        local_max_index = MIN(max_index, max_indices[type]);

        while ( ( pass ) && ( idx <= (local_max_index + lag) ) )
        {
            if ( ( pass ) && ( do_inserts ) && ( (idx + lag) >= 0 ) &&
                 ( (idx + lag) <= max_indices[type] ) &&
                 ( ((idx + lag) % 2) == 0 ) &&
                 ( ! entry_in_cache(cache_ptr, type, (idx + lag)) ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(i, %d, %d) ", type, (idx + lag));

                insert_entry(cache_ptr, type, (idx + lag), dirty_inserts,
                             H5C__NO_FLAGS_SET);
            }

            i = idx;

            while ( ( pass ) && ( i >= idx - lag ) && ( i >= 0 ) )
            {
                if ( ( pass ) && ( i >= 0 ) && ( i <= local_max_index ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(p, %d, %d) ", type, i);

                    protect_entry(cache_ptr, type, i);

                    if ( verbose )
                        HDfprintf(stdout, "(u, %d, %d) ", type, i);

                    unprotect_entry(cache_ptr, type, i, NO_CHANGE,
                                    H5C__NO_FLAGS_SET);
                }
                i--;
            }

            if ( verbose )
                HDfprintf(stdout, "\n");

            idx++;
        }
        type++;
    }

    if ( ( pass ) && ( display_stats ) ) {

        H5C_stats(cache_ptr, "test cache", display_detailed_stats);
    }

    return;

} /* hl_row_major_scan_forward() */


/*-------------------------------------------------------------------------
 * Function:	row_major_scan_backward()
 *
 * Purpose:	Do a sequence of inserts, protects, unprotects, renames,
 *		destroys while scanning backwards through the set of
 *		entries.  If pass is false on entry, do nothing.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              6/12/04
 *
 * Modifications:
 *
 * 		JRM -- 4/4/07
 * 		Added code supporting multiple read only protects.
 * 		Note that this increased the minimum lag to 10.
 *
 *-------------------------------------------------------------------------
 */

void
row_major_scan_backward(H5C_t * cache_ptr,
                        int32_t lag,
                        hbool_t verbose,
                        hbool_t reset_stats,
                        hbool_t display_stats,
                        hbool_t display_detailed_stats,
                        hbool_t do_inserts,
                        hbool_t dirty_inserts,
                        hbool_t do_renames,
                        hbool_t rename_to_main_addr,
                        hbool_t do_destroys,
			hbool_t do_mult_ro_protects,
                        int dirty_destroys,
                        int dirty_unprotects)
{
    const char * fcn_name = "row_major_scan_backward";
    int32_t type;
    int32_t idx;

    if ( verbose )
        HDfprintf(stdout, "%s(): Entering.\n", fcn_name);

    HDassert( lag >= 10 );

    type = NUMBER_OF_ENTRY_TYPES - 1;

    if ( ( pass ) && ( reset_stats ) ) {

        H5C_stats__reset(cache_ptr);
    }

    while ( ( pass ) && ( type >= 0 ) )
    {
        idx = max_indices[type] + lag;

        while ( ( pass ) && ( idx >= -lag ) )
        {
            if ( ( pass ) && ( do_inserts ) && ( (idx - lag) >= 0 ) &&
                 ( (idx - lag) <= max_indices[type] ) &&
                 ( ((idx - lag) % 2) == 1 ) &&
                 ( ! entry_in_cache(cache_ptr, type, (idx - lag)) ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(i, %d, %d) ", type, (idx - lag));

                insert_entry(cache_ptr, type, (idx - lag), dirty_inserts,
                             H5C__NO_FLAGS_SET);
            }


            if ( ( pass ) && ( (idx - lag + 1) >= 0 ) &&
                 ( (idx - lag + 1) <= max_indices[type] ) &&
                 ( ( (idx - lag + 1) % 3 ) == 0 ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(p, %d, %d) ", type, (idx - lag + 1));

                protect_entry(cache_ptr, type, (idx - lag + 1));
            }

            if ( ( pass ) && ( (idx - lag + 2) >= 0 ) &&
                 ( (idx - lag + 2) <= max_indices[type] ) &&
                 ( ( (idx - lag + 2) % 3 ) == 0 ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(u, %d, %d) ", type, (idx - lag + 2));

                unprotect_entry(cache_ptr, type, idx-lag+2, NO_CHANGE,
                                H5C__NO_FLAGS_SET);
            }


            if ( ( pass ) && ( do_renames ) && ( (idx - lag + 2) >= 0 ) &&
                 ( (idx - lag + 2) <= max_indices[type] ) &&
                 ( ( (idx - lag + 2) % 3 ) == 0 ) ) {

                rename_entry(cache_ptr, type, (idx - lag + 2),
                             rename_to_main_addr);
            }


            if ( ( pass ) && ( (idx - lag + 3) >= 0 ) &&
                 ( (idx - lag + 3) <= max_indices[type] ) &&
                 ( ( (idx - lag + 3) % 5 ) == 0 ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(p, %d, %d) ", type, (idx - lag + 3));

                protect_entry(cache_ptr, type, (idx - lag + 3));
            }

            if ( ( pass ) && ( (idx - lag + 5) >= 0 ) &&
                 ( (idx - lag + 5) <= max_indices[type] ) &&
                 ( ( (idx - lag + 5) % 5 ) == 0 ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(u, %d, %d) ", type, (idx - lag + 5));

                unprotect_entry(cache_ptr, type, idx-lag+5, NO_CHANGE,
                                H5C__NO_FLAGS_SET);
            }

	    if ( do_mult_ro_protects )
	    {
		if ( ( pass ) && ( (idx - lag + 5) >= 0 ) &&
		     ( (idx - lag + 5) < max_indices[type] ) &&
		     ( (idx - lag + 5) % 9 == 0 ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(p-ro, %d, %d) ", type,
				  (idx - lag + 5));

		    protect_entry_ro(cache_ptr, type, (idx - lag + 5));
		}

		if ( ( pass ) && ( (idx - lag + 6) >= 0 ) &&
		     ( (idx - lag + 6) < max_indices[type] ) &&
		     ( (idx - lag + 6) % 11 == 0 ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(p-ro, %d, %d) ", type,
				  (idx - lag + 6));

		    protect_entry_ro(cache_ptr, type, (idx - lag + 6));
		}

		if ( ( pass ) && ( (idx - lag + 7) >= 0 ) &&
		     ( (idx - lag + 7) < max_indices[type] ) &&
		     ( (idx - lag + 7) % 13 == 0 ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(p-ro, %d, %d) ", type,
				  (idx - lag + 7));

		    protect_entry_ro(cache_ptr, type, (idx - lag + 7));
		}

		if ( ( pass ) && ( (idx - lag + 7) >= 0 ) &&
		     ( (idx - lag + 7) < max_indices[type] ) &&
		     ( (idx - lag + 7) % 9 == 0 ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(u-ro, %d, %d) ", type,
				  (idx - lag + 7));

		    unprotect_entry(cache_ptr, type, (idx - lag + 7),
				    FALSE, H5C__NO_FLAGS_SET);
		}

		if ( ( pass ) && ( (idx - lag + 8) >= 0 ) &&
		     ( (idx - lag + 8) < max_indices[type] ) &&
		     ( (idx - lag + 8) % 11 == 0 ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(u-ro, %d, %d) ", type,
				  (idx - lag + 8));

		    unprotect_entry(cache_ptr, type, (idx - lag + 8),
				    FALSE, H5C__NO_FLAGS_SET);
		}

		if ( ( pass ) && ( (idx - lag + 9) >= 0 ) &&
		     ( (idx - lag + 9) < max_indices[type] ) &&
		     ( (idx - lag + 9) % 13 == 0 ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(u-ro, %d, %d) ", type,
				  (idx - lag + 9));

		    unprotect_entry(cache_ptr, type, (idx - lag + 9),
				    FALSE, H5C__NO_FLAGS_SET);
		}
	    } /* if ( do_mult_ro_protects ) */

            if ( ( pass ) && ( idx >= 0 ) && ( idx <= max_indices[type] ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(p, %d, %d) ", type, idx);

                protect_entry(cache_ptr, type, idx);
            }


            if ( ( pass ) && ( (idx + lag - 2) >= 0 ) &&
                 ( (idx + lag - 2) <= max_indices[type] ) &&
                 ( ( (idx + lag - 2) % 7 ) == 0 ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(u, %d, %d) ", type, (idx + lag - 2));

                unprotect_entry(cache_ptr, type, idx+lag-2, NO_CHANGE,
                                H5C__NO_FLAGS_SET);
            }

            if ( ( pass ) && ( (idx + lag - 1) >= 0 ) &&
                 ( (idx + lag - 1) <= max_indices[type] ) &&
                 ( ( (idx + lag - 1) % 7 ) == 0 ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(p, %d, %d) ", type, (idx + lag - 1));

                protect_entry(cache_ptr, type, (idx + lag - 1));
            }


            if ( do_destroys ) {

                if ( ( pass ) && ( (idx + lag) >= 0 ) &&
                     ( ( idx + lag) <= max_indices[type] ) ) {

                    switch ( (idx + lag) %4 ) {

                        case 0:
                            if ( (entries[type])[idx+lag].is_dirty ) {

                                unprotect_entry(cache_ptr, type, idx + lag,
                                                NO_CHANGE, H5C__NO_FLAGS_SET);
                            } else {

                                unprotect_entry(cache_ptr, type, idx + lag,
                                                dirty_unprotects,
                                                H5C__NO_FLAGS_SET);
                            }
                            break;

                        case 1: /* we just did an insert */
                            unprotect_entry(cache_ptr, type, idx + lag,
                                            NO_CHANGE, H5C__NO_FLAGS_SET);
                            break;

                        case 2:
                            if ( (entries[type])[idx + lag].is_dirty ) {

                                unprotect_entry(cache_ptr, type, idx + lag,
                                                NO_CHANGE, H5C__DELETED_FLAG);
                            } else {

                                unprotect_entry(cache_ptr, type, idx + lag,
                                                dirty_destroys,
                                                H5C__DELETED_FLAG);
                            }
                            break;

                        case 3: /* we just did an insrt */
                            unprotect_entry(cache_ptr, type, idx + lag,
                                            NO_CHANGE, H5C__DELETED_FLAG);
                            break;

                        default:
                            HDassert(0); /* this can't happen... */
                            break;
                    }
                }
            } else {

                if ( ( pass ) && ( (idx + lag) >= 0 ) &&
                     ( ( idx + lag) <= max_indices[type] ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(u, %d, %d) ", type, (idx - lag));

                    unprotect_entry(cache_ptr, type, idx + lag,
                                    dirty_unprotects, H5C__NO_FLAGS_SET);
                }
            }

            if ( verbose )
                HDfprintf(stdout, "\n");

            idx--;
        }
        type--;
    }

    if ( ( pass ) && ( display_stats ) ) {

        H5C_stats(cache_ptr, "test cache", display_detailed_stats);
    }

    return;

} /* row_major_scan_backward() */


/*-------------------------------------------------------------------------
 * Function:	hl_row_major_scan_backward()
 *
 * Purpose:	Do a high locality sequence of inserts, protects, and
 *		unprotects while scanning through the set of entries.
 *		If pass is false on entry, do nothing.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              10/21/04
 *
 * Modifications:
 *
 *		JRM -- 1/21/05
 *		Added the max_index parameter to allow the caller to
 *		throttle the size of the inner loop, and thereby the
 *		execution time of the function.
 *
 *-------------------------------------------------------------------------
 */

void
hl_row_major_scan_backward(H5C_t * cache_ptr,
                           int32_t max_index,
                           hbool_t verbose,
                           hbool_t reset_stats,
                           hbool_t display_stats,
                           hbool_t display_detailed_stats,
                           hbool_t do_inserts,
                           hbool_t dirty_inserts)
{
    const char * fcn_name = "hl_row_major_scan_backward";
    int32_t type;
    int32_t idx;
    int32_t i;
    int32_t lag = 100;
    int32_t local_max_index;

    if ( verbose )
        HDfprintf(stdout, "%s(): entering.\n", fcn_name);

    HDassert( lag > 5 );
    HDassert( max_index >= 200 );
    HDassert( max_index <= MAX_ENTRIES );

    type = NUMBER_OF_ENTRY_TYPES - 1;

    if ( ( pass ) && ( reset_stats ) ) {

        H5C_stats__reset(cache_ptr);
    }

    while ( ( pass ) && ( type >= 0 ) )
    {
        idx = max_indices[type] + lag;

        local_max_index = MIN(max_index, max_indices[type]);

        while ( ( pass ) && ( idx >= -lag ) )
        {
            if ( ( pass ) && ( do_inserts ) && ( (idx + lag) >= 0 ) &&
                 ( (idx + lag) <= local_max_index ) &&
                 ( ((idx + lag) % 2) == 0 ) &&
                 ( ! entry_in_cache(cache_ptr, type, (idx + lag)) ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(i, %d, %d) ", type, (idx + lag));

                insert_entry(cache_ptr, type, (idx + lag), dirty_inserts,
                             H5C__NO_FLAGS_SET);
            }

            i = idx;

            while ( ( pass ) && ( i >= idx - lag ) && ( i >= 0 ) )
            {
                if ( ( pass ) && ( i >= 0 ) && ( i <= local_max_index ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(p, %d, %d) ", type, i);

                    protect_entry(cache_ptr, type, i);

                    if ( verbose )
                        HDfprintf(stdout, "(u, %d, %d) ", type, i);

                    unprotect_entry(cache_ptr, type, i, NO_CHANGE,
                                    H5C__NO_FLAGS_SET);
                }
                i--;
            }

            if ( verbose )
                HDfprintf(stdout, "\n");

            idx--;
        }
        type--;
    }

    if ( ( pass ) && ( display_stats ) ) {

        H5C_stats(cache_ptr, "test cache", display_detailed_stats);
    }

    return;

} /* hl_row_major_scan_backward() */


/*-------------------------------------------------------------------------
 * Function:	col_major_scan_forward()
 *
 * Purpose:	Do a sequence of inserts, protects, and unprotects
 *		while scanning through the set of entries.  If
 *		pass is false on entry, do nothing.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              6/23/04
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */

void
col_major_scan_forward(H5C_t * cache_ptr,
                       int32_t lag,
                       hbool_t verbose,
                       hbool_t reset_stats,
                       hbool_t display_stats,
                       hbool_t display_detailed_stats,
                       hbool_t do_inserts,
                       hbool_t dirty_inserts,
                       int dirty_unprotects)
{
    const char * fcn_name = "col_major_scan_forward()";
    int32_t type;
    int32_t idx;

    if ( verbose )
        HDfprintf(stdout, "%s: entering.\n", fcn_name);

    HDassert( lag > 5 );

    type = 0;

    if ( ( pass ) && ( reset_stats ) ) {

        H5C_stats__reset(cache_ptr);
    }

    idx = -lag;

    while ( ( pass ) && ( (idx - lag) <= MAX_ENTRIES ) )
    {
        type = 0;

        while ( ( pass ) && ( type < NUMBER_OF_ENTRY_TYPES ) )
        {
            if ( ( pass ) && ( do_inserts ) && ( (idx + lag) >= 0 ) &&
                 ( (idx + lag) <= max_indices[type] ) &&
                 ( ((idx + lag) % 3) == 0 ) &&
                 ( ! entry_in_cache(cache_ptr, type, (idx + lag)) ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(i, %d, %d) ", type, (idx + lag));

                insert_entry(cache_ptr, type, (idx + lag), dirty_inserts,
                             H5C__NO_FLAGS_SET);
            }

            if ( ( pass ) && ( idx >= 0 ) && ( idx <= max_indices[type] ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(p, %d, %d) ", type, idx);

                protect_entry(cache_ptr, type, idx);
            }

            if ( ( pass ) && ( (idx - lag) >= 0 ) &&
                 ( (idx - lag) <= max_indices[type] ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(u, %d, %d) ", type, (idx - lag));

                unprotect_entry(cache_ptr, type, idx - lag,
                                dirty_unprotects, H5C__NO_FLAGS_SET);
            }

            if ( verbose )
                HDfprintf(stdout, "\n");

            type++;
        }

        idx++;
    }

    if ( ( pass ) && ( display_stats ) ) {

        H5C_stats(cache_ptr, "test cache", display_detailed_stats);
    }

    return;

} /* col_major_scan_forward() */


/*-------------------------------------------------------------------------
 * Function:	hl_col_major_scan_forward()
 *
 * Purpose:	Do a high locality sequence of inserts, protects, and
 *		unprotects while scanning through the set of entries.  If
 *		pass is false on entry, do nothing.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              19/25/04
 *
 * Modifications:
 *
 *		JRM -- 1/21/05
 *		Added the max_index parameter to allow the caller to
 *		throttle the size of the inner loop, and thereby the
 *		execution time of the function.
 *
 *-------------------------------------------------------------------------
 */

void
hl_col_major_scan_forward(H5C_t * cache_ptr,
                          int32_t max_index,
                          hbool_t verbose,
                          hbool_t reset_stats,
                          hbool_t display_stats,
                          hbool_t display_detailed_stats,
                          hbool_t do_inserts,
                          hbool_t dirty_inserts,
                          int dirty_unprotects)
{
    const char * fcn_name = "hl_col_major_scan_forward()";
    int32_t type;
    int32_t idx;
    int32_t lag = 200;
    int32_t i;
    int32_t local_max_index;

    if ( verbose )
        HDfprintf(stdout, "%s: entering.\n", fcn_name);

    HDassert( lag > 5 );
    HDassert( max_index >= 500 );
    HDassert( max_index <= MAX_ENTRIES );

    type = 0;

    if ( ( pass ) && ( reset_stats ) ) {

        H5C_stats__reset(cache_ptr);
    }

    idx = 0;

    local_max_index = MIN(max_index, MAX_ENTRIES);

    while ( ( pass ) && ( idx <= local_max_index ) )
    {

        i = idx;

        while ( ( pass ) && ( i >= 0 ) && ( i >= (idx - lag) ) ) {

            type = 0;

            while ( ( pass ) && ( type < NUMBER_OF_ENTRY_TYPES ) )
            {
                if ( ( pass ) && ( do_inserts ) && ( i == idx ) &&
                     ( i <= local_max_index ) &&
                     ( (i % 3) == 0 ) &&
                     ( ! entry_in_cache(cache_ptr, type, i) ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(i, %d, %d) ", type, i);

                    insert_entry(cache_ptr, type, i, dirty_inserts,
                                 H5C__NO_FLAGS_SET);
                }

                if ( ( pass ) && ( i >= 0 ) && ( i <= local_max_index ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(p, %d, %d) ", type, i);

                    protect_entry(cache_ptr, type, i);
                }

                if ( ( pass ) && ( i >= 0 ) &&
                     ( i <= max_indices[type] ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(u, %d, %d) ", type, i);

                    unprotect_entry(cache_ptr, type, i,
                                    dirty_unprotects, H5C__NO_FLAGS_SET);
                }

                if ( verbose )
                    HDfprintf(stdout, "\n");

                type++;
            }

            i--;
        }

        idx++;
    }

    if ( ( pass ) && ( display_stats ) ) {

        H5C_stats(cache_ptr, "test cache", display_detailed_stats);
    }

    return;

} /* hl_col_major_scan_forward() */


/*-------------------------------------------------------------------------
 * Function:	col_major_scan_backward()
 *
 * Purpose:	Do a sequence of inserts, protects, and unprotects
 *		while scanning backwards through the set of
 *		entries.  If pass is false on entry, do nothing.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              6/23/04
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */

void
col_major_scan_backward(H5C_t * cache_ptr,
                        int32_t lag,
                        hbool_t verbose,
                        hbool_t reset_stats,
                        hbool_t display_stats,
                        hbool_t display_detailed_stats,
                        hbool_t do_inserts,
                        hbool_t dirty_inserts,
                        int dirty_unprotects)
{
    const char * fcn_name = "col_major_scan_backward()";
    int mile_stone = 1;
    int32_t type;
    int32_t idx;

    if ( verbose )
        HDfprintf(stdout, "%s: entering.\n", fcn_name);

    HDassert( lag > 5 );

    if ( ( pass ) && ( reset_stats ) ) {

        H5C_stats__reset(cache_ptr);
    }

    idx = MAX_ENTRIES + lag;

    if ( verbose ) /* 1 */
        HDfprintf(stdout, "%s: point %d.\n", fcn_name, mile_stone++);


    while ( ( pass ) && ( (idx + lag) >= 0 ) )
    {
        type = NUMBER_OF_ENTRY_TYPES - 1;

        while ( ( pass ) && ( type >= 0 ) )
        {
            if ( ( pass ) && ( do_inserts) && ( (idx - lag) >= 0 ) &&
                 ( (idx - lag) <= max_indices[type] ) &&
                 ( ((idx - lag) % 3) == 0 ) &&
                 ( ! entry_in_cache(cache_ptr, type, (idx - lag)) ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(i, %d, %d) ", type, (idx - lag));

                insert_entry(cache_ptr, type, (idx - lag), dirty_inserts,
                             H5C__NO_FLAGS_SET);
            }

            if ( ( pass ) && ( idx >= 0 ) && ( idx <= max_indices[type] ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(p, %d, %d) ", type, idx);

                protect_entry(cache_ptr, type, idx);
            }

            if ( ( pass ) && ( (idx + lag) >= 0 ) &&
                 ( (idx + lag) <= max_indices[type] ) ) {

                if ( verbose )
                    HDfprintf(stdout, "(u, %d, %d) ", type, (idx + lag));

                unprotect_entry(cache_ptr, type, idx + lag,
                                dirty_unprotects, H5C__NO_FLAGS_SET);
            }

            if ( verbose )
                HDfprintf(stdout, "\n");

            type--;
        }

        idx--;
    }

    if ( verbose ) /* 2 */
        HDfprintf(stdout, "%s: point %d.\n", fcn_name, mile_stone++);

    if ( ( pass ) && ( display_stats ) ) {

        H5C_stats(cache_ptr, "test cache", display_detailed_stats);
    }

    if ( verbose )
        HDfprintf(stdout, "%s: exiting.\n", fcn_name);

    return;

} /* col_major_scan_backward() */


/*-------------------------------------------------------------------------
 * Function:	hl_col_major_scan_backward()
 *
 * Purpose:	Do a high locality sequence of inserts, protects, and
 *		unprotects while scanning backwards through the set of
 *		entries.  If pass is false on entry, do nothing.
 *
 * Return:	void
 *
 * Programmer:	John Mainzer
 *              10/25/04
 *
 * Modifications:
 *
 *		JRM -- 1/21/05
 *		Added the max_index parameter to allow the caller to
 *		throttle the size of the inner loop, and thereby the
 *		execution time of the function.
 *
 *-------------------------------------------------------------------------
 */

void
hl_col_major_scan_backward(H5C_t * cache_ptr,
                           int32_t max_index,
                           hbool_t verbose,
                           hbool_t reset_stats,
                           hbool_t display_stats,
                           hbool_t display_detailed_stats,
                           hbool_t do_inserts,
                           hbool_t dirty_inserts,
                           int dirty_unprotects)
{
    const char * fcn_name = "hl_col_major_scan_backward()";
    int32_t type;
    int32_t idx;
    int32_t lag = 50;
    int32_t i;
    int32_t local_max_index;

    if ( verbose )
        HDfprintf(stdout, "%s: entering.\n", fcn_name);

    HDassert( lag > 5 );
    HDassert( max_index >= 500 );
    HDassert( max_index <= MAX_ENTRIES );

    type = 0;

    local_max_index = MIN(max_index, MAX_ENTRIES);

    if ( ( pass ) && ( reset_stats ) ) {

        H5C_stats__reset(cache_ptr);
    }

    idx = local_max_index;

    while ( ( pass ) && ( idx >= 0 ) )
    {

        i = idx;

        while ( ( pass ) && ( i <= local_max_index ) && ( i <= (idx + lag) ) ) {

            type = 0;

            while ( ( pass ) && ( type < NUMBER_OF_ENTRY_TYPES ) )
            {
                if ( ( pass ) && ( do_inserts ) && ( i == idx ) &&
                     ( i <= local_max_index ) &&
                     ( ! entry_in_cache(cache_ptr, type, i) ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(i, %d, %d) ", type, i);

                    insert_entry(cache_ptr, type, i, dirty_inserts,
                                 H5C__NO_FLAGS_SET);
                }

                if ( ( pass ) && ( i >= 0 ) && ( i <= local_max_index ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(p, %d, %d) ", type, i);

                    protect_entry(cache_ptr, type, i);
                }

                if ( ( pass ) && ( i >= 0 ) &&
                     ( i <= local_max_index ) ) {

                    if ( verbose )
                        HDfprintf(stdout, "(u, %d, %d) ", type, i);

                    unprotect_entry(cache_ptr, type, i,
                                    dirty_unprotects, H5C__NO_FLAGS_SET);
                }

                if ( verbose )
                    HDfprintf(stdout, "\n");

                type++;
            }

            i++;
        }

        idx--;
    }

    if ( ( pass ) && ( display_stats ) ) {

        H5C_stats(cache_ptr, "test cache", display_detailed_stats);
    }

    return;

} /* hl_col_major_scan_backward() */


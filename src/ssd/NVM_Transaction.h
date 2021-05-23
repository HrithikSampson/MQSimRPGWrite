#ifndef NVM_TRANSACTION_H
#define NVM_TRANSACTION_H

#include<iostream>
#include <list>
#include "../sim/Sim_Defs.h"
#include "../sim/Engine.h"
#include "User_Request.h"
//#define Update_transaction_statistics( transaction )  Up_trx_stats_real( __FUNCTION__ , __FILE__ , __LINE__ , transaction )
namespace SSD_Components
{
	class User_Request;
	
	enum class Transaction_Type { READ, WRITE, ERASE, UNKOWN };
	enum class Transaction_Source_Type { USERIO, CACHE, GC_WL, MAPPING };

	class NVM_Transaction
	{
	public:
		NVM_Transaction(stream_id_type stream_id, Transaction_Source_Type source, Transaction_Type type, User_Request* user_request) :
			Stream_id(stream_id), Source(source), Type(type), UserIORequest(user_request), Issue_time(Simulator->Time()), STAT_transfer_time(TRANSFER_TIME) {
            Stream_id=0;
			if(type==Transaction_Type::READ)
			{
				STAT_execution_time=READ_TIME;
			}
			else if(type==Transaction_Type::WRITE)
			{
				STAT_execution_time=WRITE_TIME;
				
			}
			else
			{
				STAT_execution_time=INVALID_TIME/5;	
				
			}
            RPG_WR=false;

		}
		bool RPG_WR;
		stream_id_type Stream_id;
		Transaction_Source_Type Source;
		Transaction_Type Type;
		User_Request* UserIORequest;
		//std::list<NVM_Transaction*>::iterator RelatedNodeInQueue;//Just used for high performance linkedlist insertion/deletion

		sim_time_type Issue_time;
		/* Used to calculate service time and transfer time for a normal read/program operation used to respond to the host IORequests.
		In other words, these variables are not important if FlashTransactions is used for garbage collection.*/
		sim_time_type STAT_execution_time, STAT_transfer_time;
	};
}

#endif //!NVM_TRANSACTION_H

#ifndef __SCAN_CONTEXT_H
#define __SCAN_CONTEXT_H

#include <iostream>
#include <Eigen/Dense>

#include "data_type.h"

#define NUM_RING 20
#define NUM_SECTOR 60

namespace context 
{


template<typename T>
static constexpr T MAX_RADIUS()
{
	return static_cast<T>(10.0);
}
	
template<typename T>
class ScanContext
{
public:
	using DataType = T;

	ScanContext()
	{

	}

	const Eigen::Matrix<DataType, NUM_RING, NUM_SECTOR>& getDesc() const
        {
                return desc;
        }
	
	const Eigen::Matrix<DataType, NUM_RING, 1> caculateResult( const sensor::LaserScan &scan )
        {
                makeScanContext( scan );
	
		Eigen::Matrix<DataType, NUM_RING, 1> mean = desc.rowwise().mean();
		return mean;
	}
	
	void makeScanContext( const sensor::LaserScan &scan )
        {
                int ring_idx = 0, sctor_idx = 0;
                DataType radians = -3.141592653;

                for( int i = 0; i < scan.size(); i ++ ){
                        DataType dist = scan.ranges[i];
                        DataType angle = rad2deg<DataType>( radians ) + 180.0f;

                        if( dist >= 0.009999998 && dist <= MAX_RADIUS<DataType>() ){
                                ring_idx = std::max( std::min( NUM_RING - 1, static_cast<int>( ceil( ( dist / MAX_RADIUS<DataType>() ) * NUM_RING ) ) ), 0 );

                                sctor_idx = std::max( std::min( NUM_SECTOR - 1, static_cast<int>(ceil( ( angle / 360.0) * NUM_SECTOR ) ) ), 0 );


                                desc( ring_idx, sctor_idx ) += 1;
                        }

                        radians += 0.0043633231;
                }
        }

private:
	template<typename TT>
        const TT rad2deg( const TT radians)
        {
                return radians * 180.0 / M_PI;
        }

        template<typename TT>
        const TT deg2rad( const TT angle )
        {
                return angle * M_PI / 180.0;
        }


private:

	Eigen::Matrix<DataType, NUM_RING, NUM_SECTOR> desc = Eigen::Matrix<DataType, NUM_RING, NUM_SECTOR>::Zero();
	Eigen::Matrix<DataType, NUM_SECTOR, NUM_SECTOR> covarince = Eigen::Matrix<DataType, NUM_SECTOR, NUM_SECTOR>::Zero();
};

}

#endif

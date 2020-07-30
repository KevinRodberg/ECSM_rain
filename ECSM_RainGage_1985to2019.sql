/* Formatted on 4/13/2020 1:21:48 PM (QP5 v5.287) */
  SELECT STATION,AGENCY,XCOORD,YCOORD,daily.DBKEY,
         daily.daily_date,
         daily.VALUE,
         daily.code
    FROM dmdbase.dm_daily_data daily,
         (  SELECT *
              FROM keyword_tab k
             WHERE     data_type = 'RAIN'
                   AND FREQUENCY = 'DA'
                   AND end_date > TO_DATE ('01-01-1984', 'MM-DD-YYYY')
                   AND k.xcoord > 672000
                   AND k.xcoord < 985000
                   AND k.ycoord > 143000
                   AND k.ycoord < 1202000
          ORDER BY xcoord, ycoord) key
   WHERE     key.dbkey = daily.dbkey
         AND daily.daily_date > TO_DATE ('12-31-1984', 'MM-DD-YYYY')
         AND daily.daily_date <= TO_DATE ('12-31-2019', 'MM-DD-YYYY')
ORDER BY key.station, key.dbkey, daily.daily_date
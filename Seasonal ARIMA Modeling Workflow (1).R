library(fpp3)

# Read Data and Make it a tsibble

X <- read.csv("US Electricity.csv") %>%
  transmute(DATE = yearmonth(DATE),
            ELEC = ELEC) %>%
  as_tsibble(index = DATE) %>%
  filter(DATE >= yearmonth("2010 Jan"))

X %>% autoplot()

# Define Training and Testing sets

TR <- X %>% filter(DATE <= yearmonth("2019 Aug"))
TE <- X %>% filter(DATE >= yearmonth("2019 Sep"))


# Test stationarity and 
# "difference" training set

X %>% 
  mutate(d1 = difference(ELEC),
         d12 = difference(ELEC,12),
         d1.12 = difference(d1,12)) -> TR

TR %>% ACF(ELEC) %>% autoplot()
TR %>% ACF(d1) %>% autoplot()
TR %>% ACF(d12) %>% autoplot()
TR %>% ACF(d1.12) %>% autoplot()


# ADF test Ho: Data is non-stationary (Want low p-value)
TR$ELEC %>% tseries::adf.test()
# KPSS H0: Trend in data is deterministic (We want a very high p-value)
TR %>% features(ELEC, unitroot_kpss)
TR$d1 %>% na.omit() %>% tseries::adf.test()
TR %>% features(d1, unitroot_kpss)
TR$d12 %>% na.omit() %>% tseries::adf.test()
TR %>% features(d12, unitroot_kpss)
TR$d1.12 %>% na.omit() %>% tseries::adf.test()
TR %>% features(d1.12, unitroot_kpss)

X %>% features(ELEC, unitroot_ndiffs)


# d=1, D=0 p, q, P, Q
# d=1, D=1

# Examine ACF and PACF to guess possible models


TR %>% gg_tsdisplay(ELEC,"partial", lag_max = 48)
# q = Q = 0, p = 5, P = 0

TR %>% gg_tsdisplay(d1,"partial", lag_max = 48)
# q = Q = 0, p = 4, P = 0

TR %>% gg_tsdisplay(d12,"partial", lag_max = 48)
# q = Q = 0, p = 1, P = 0

TR %>% gg_tsdisplay(d1.12,"partial", lag_max = 48)
# p = P = 0, q = 1, Q = 0

# Use Trial and Error to fit several models

m <- TR %>%
  model(mg1 = ARIMA(ELEC ~  pdq(0,1,1) + PDQ(0,1,0)),
        mg2 = ARIMA(ELEC ~  pdq(1,0,0) + PDQ(0,1,0)),
        ma = ARIMA(ELEC),
        mets = ETS(ELEC))

m %>% glance()

m %>% select(mg1) %>% report()
m %>% select(mg2) %>% report()
m %>% select(ma) %>% report()
m %>% select(mets) %>% report()


# Validate residual independence assumption

m %>% augment() %>%
  features(.resid, ljung_box, lag = 24)

m %>% select(ma) %>% gg_tsresiduals() 
m %>% select(mg1) %>% gg_tsresiduals()
m %>% select(mg2) %>% gg_tsresiduals()
m %>% select(mets) %>% gg_tsresiduals()

# Examine accuracy metrics

f <- m %>% forecast(h=24)

rbind(m %>% accuracy(),
      f %>% accuracy(TE))
      


# Cross Validation
# Read Section 5.10 of FPP3
#
X.CV <- X %>%
  filter(DATE >= yearmonth("2010 Jan")) %>%
  filter(DATE <= yearmonth("2019 Aug")) %>% 
  stretch_tsibble(.init = 72, .step = 1)


mC <- X.CV %>% 
  model(mg1 = ARIMA(ELEC ~  pdq(0,1,1) + PDQ(0,1,1)),
        mg2 = ARIMA(ELEC ~ pdq(1,0,0) + PDQ(0,1,0)),
        ma =  ARIMA(ELEC ~ pdq(1,0,0) + PDQ(2,1,0)),
        mets = ETS(ELEC ~ error("M") + trend("N") + season("A")))




mC %>%
  forecast(h = 24) %>%
  group_by(.id, .model) %>%
  mutate(h = row_number()) %>%
  ungroup %>% 
  as_fable(response = "ELEC", distribution = ELEC)-> fCV

fCV %>%
  accuracy(X, by = c("h", ".model")) %>%
  ggplot(aes(x = h, y = MAPE, color = .model)) +
  geom_line()


fCV %>% filter(.id == 45, .model == "ma") %>% 
  autoplot() + 
  geom_line(aes(y = ELEC), data = X)

m %>% select(ma) %>% forecast(h=24) %>% autoplot(TR)+
  geom_point(aes(y = ELEC), data = TE, col = "red")




# Correlation of inter-model residuals

m %>% 
  augment() %>%
  as_tibble() %>%
  select(DATE,.model,.resid) %>%
  spread(.model,.resid) %>%
  select(-DATE) %>%
  cor()
